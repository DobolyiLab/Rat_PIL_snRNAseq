# =============================================================================
# 00_setup.R  -- packages, options, shared helper functions
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(SoupX)
  library(scDblFinder); library(SingleCellExperiment)
  library(Matrix); library(data.table); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(future)
  library(readr); library(yaml)
})

set.seed(42)
plan("sequential")
options(future.globals.maxSize = 20 * 1024^3, future.seed = TRUE)

params <- yaml::read_yaml("config/params.yml")

# S4Vectors (via SingleCellExperiment) and AnnotationDbi export functions with
# the same names as dplyr verbs. Depending on attach order they can mask dplyr,
# and the failure is confusing: S4Vectors::rename() evaluates its arguments and
# reports "object 'x' not found". Pin the dplyr versions explicitly.
rename    <- dplyr::rename
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
arrange   <- dplyr::arrange
count     <- dplyr::count
first     <- dplyr::first
slice     <- dplyr::slice

for (d in c("figures", "data/processed", "data/processed/qc")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
source("utils/plotting.R")

SAMPLE_IDS <- vapply(params$samples, function(s) s$id, character(1))

# -----------------------------------------------------------------------------
# Matrix input
# The data are NOT in cellranger `outs/` layout. They sit in one flat folder
# with per-sample prefixes (R1_Cont_matrix.mtx, R1_Cont_barcodes.tsv, ...),
# so Read10X() cannot be used -- ReadMtx() with explicit file paths is correct.
# -----------------------------------------------------------------------------
.pick_file <- function(files, prefix, pattern) {
  hit <- files[grepl(paste0("^", prefix), files) & grepl(pattern, files, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_matrix_triplet <- function(dir, prefix) {
  fs <- list.files(dir)
  tri <- list(
    mtx      = .pick_file(fs, prefix, "matrix\\.mtx(\\.gz)?$"),
    barcodes = .pick_file(fs, prefix, "barcodes\\.tsv(\\.gz)?$"),
    features = .pick_file(fs, prefix, "(features|genes)\\.tsv(\\.gz)?$")
  )
  # Fallback: standard cellranger layout, where the files carry no sample
  # prefix (this is how the raw_feature_bc_matrix folders are organized).
  if (anyNA(unlist(tri))) {
    tri <- list(
      mtx      = .pick_file(fs, "", "^matrix\\.mtx(\\.gz)?$"),
      barcodes = .pick_file(fs, "", "^barcodes\\.tsv(\\.gz)?$"),
      features = .pick_file(fs, "", "^(features|genes)\\.tsv(\\.gz)?$")
    )
  }
  tri
}

.n_columns <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  on.exit(close(con))
  length(strsplit(readLines(con, n = 1), "\t")[[1]])
}

# Reads one sample. `which` is "filtered" or "raw"; the raw folder is optional
# and only used by SoupX in step 03.
read_sample_matrix <- function(sample, which = "filtered") {
  dir <- if (which == "filtered") params$matrix_dir else sample$raw_dir
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) {
    stop(sprintf("%s matrix directory for %s not available (%s)",
                 which, sample$id, if (is.null(dir)) "not set in params.yml" else dir),
         call. = FALSE)
  }
  prefix <- if (!is.null(sample$prefix)) sample$prefix else sample$id
  tri <- find_matrix_triplet(dir, prefix)
  
  if (anyNA(unlist(tri))) {
    stop(sprintf(
      "Incomplete matrix triplet for '%s' in %s\n  matrix  : %s\n  barcodes: %s\n  features: %s\n  files present: %s",
      prefix, dir, tri$mtx, tri$barcodes, tri$features,
      paste(list.files(dir), collapse = ", ")), call. = FALSE)
  }
  
  fp <- file.path(dir, tri$features)
  feat_col <- if (.n_columns(fp) >= 2) 2 else 1
  
  m <- ReadMtx(mtx = file.path(dir, tri$mtx),
               cells = file.path(dir, tri$barcodes),
               features = fp,
               feature.column = feat_col,
               strip.suffix = FALSE)
  
  message(sprintf("[read] %s (%s): %d features x %d cells (feature.column = %d)",
                  sample$id, which, nrow(m), ncol(m), feat_col))
  m
}

# TRUE only if every sample has a usable raw droplet matrix (SoupX prerequisite).
raw_matrices_available <- function() {
  all(vapply(params$samples, function(s) {
    d <- s$raw_dir
    if (is.null(d) || !nzchar(d) || !dir.exists(d)) return(FALSE)
    prefix <- if (!is.null(s$prefix)) s$prefix else s$id
    !anyNA(unlist(find_matrix_triplet(d, prefix)))
  }, logical(1)))
}

# Fail once, up front, with every problem named.
check_sample_paths <- function(need_raw = FALSE) {
  if (is.null(params$matrix_dir) || !dir.exists(params$matrix_dir)) {
    stop("params$matrix_dir does not exist: ",
         if (is.null(params$matrix_dir)) "<not set>" else params$matrix_dir,
         "\n(working directory: ", getwd(), ")", call. = FALSE)
  }
  bad <- character(0)
  for (s in params$samples) {
    prefix <- if (!is.null(s$prefix)) s$prefix else s$id
    tri <- find_matrix_triplet(params$matrix_dir, prefix)
    if (anyNA(unlist(tri))) {
      bad <- c(bad, sprintf("  %s: missing %s", s$id,
                            paste(names(tri)[is.na(unlist(tri))], collapse = ", ")))
    }
  }
  if (length(bad)) {
    stop("Incomplete filtered matrices in ", params$matrix_dir, ":\n",
         paste(bad, collapse = "\n"),
         "\nFiles present: ", paste(list.files(params$matrix_dir), collapse = ", "),
         call. = FALSE)
  }
  if (need_raw && !raw_matrices_available()) {
    warning("Raw (unfiltered) droplet matrices are not available for all samples. ",
            "SoupX cannot be run; step 03 will fall back to ",
            params$ambient$method_without_raw, ". See README.md.")
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Seurat v4 / v5 compatibility for the counts/data slots
# -----------------------------------------------------------------------------
.seurat_v5 <- utils::packageVersion("SeuratObject") >= "5.0.0"

# In Seurat v5 an assay can hold one layer PER SAMPLE (counts.R1_Cont, ...).
# LayerData() then silently returns only the first one -- which is how a
# 29,325-cell object turns into an 8,246-column matrix. Join them explicitly.
split_layers <- function(obj, assay = "RNA", slot = "counts") {
  if (!.seurat_v5) return(character(0))
  lyr <- SeuratObject::Layers(obj[[assay]])
  grep(paste0("^", slot, "\\."), lyr, value = TRUE)
}

join_layers_if_split <- function(obj, assay = "RNA") {
  if (!.seurat_v5) return(obj)
  if (length(split_layers(obj, assay, "counts")) == 0 &&
      length(split_layers(obj, assay, "data")) == 0) return(obj)
  message("[layers] joining split ", assay, " layers: ",
          paste(SeuratObject::Layers(obj[[assay]]), collapse = ", "))
  obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
  obj
}

get_expr <- function(obj, assay = "RNA", slot = "data") {
  if (!.seurat_v5) return(Seurat::GetAssayData(obj, assay = assay, slot = slot))
  
  spl <- split_layers(obj, assay, slot)
  if (length(spl) > 1) {
    stop(sprintf(paste0("Assay '%s' has %d split layers for '%s' (%s). LayerData() would ",
                        "return only the first one and silently drop most cells.\n",
                        "Call obj <- join_layers_if_split(obj, \"%s\") first."),
                 assay, length(spl), slot, paste(spl, collapse = ", "), assay),
         call. = FALSE)
  }
  m <- SeuratObject::LayerData(obj, assay = assay, layer = slot)
  if (ncol(m) != ncol(obj)) {
    stop(sprintf("Extracted matrix has %d columns but the object has %d cells (assay '%s', layer '%s').",
                 ncol(m), ncol(obj), assay, slot), call. = FALSE)
  }
  m
}

# TRUE if the RNA `data` layer does not hold log-normalized values yet.
# Two distinct situations both mean "not normalized":
#   * Seurat v5 leaves `data` EMPTY after CreateSeuratObject() until
#     NormalizeData() runs -- it does not mirror counts the way v4 did
#   * `data` exists but is byte-identical to `counts`
# The original `max(data) <= 0` test caught neither.
rna_is_unnormalized <- function(obj, assay = "RNA") {
  d <- tryCatch({
    if (.seurat_v5) {
      if (!any(grepl("^data", SeuratObject::Layers(obj[[assay]])))) return(TRUE)
      SeuratObject::LayerData(obj, assay = assay, layer = "data")
    } else {
      Seurat::GetAssayData(obj, assay = assay, slot = "data")
    }
  }, error = function(e) NULL)
  
  if (is.null(d) || nrow(d) == 0 || ncol(d) == 0) return(TRUE)   # empty layer
  
  cts <- get_expr(obj, assay, "counts")
  if (!identical(dim(d), dim(cts))) return(FALSE)
  n <- min(2000, length(d@x))
  if (n == 0) return(TRUE)
  idx <- seq_len(n)
  isTRUE(all.equal(as.numeric(d@x[idx]), as.numeric(cts@x[idx])))
}

ensure_rna_normalized <- function(obj, assay = "RNA") {
  obj <- join_layers_if_split(obj, assay)   # safe to call even if already joined
  if (rna_is_unnormalized(obj, assay)) {
    message("RNA assay is not log-normalized -> running NormalizeData()")
    obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
    obj <- join_layers_if_split(obj, assay) # v5 may re-split `data` per sample
  }
  obj
}

# -----------------------------------------------------------------------------
# Mitochondrial genes -- RAT, and the naming depends on the annotation source:
#   Ensembl (Rnor_6.0 / mRatBN7.2) : "Mt-nd1", "Mt-co1"
#   10x prebuilt rat reference     : "AY172581.x" (accession-based)
#   NCBI RefSeq (mRatBN7.2)        : bare "ND1", "COX1", "CYTB"  <-- this dataset
# The RefSeq names are anchored with ^...$ so nuclear genes such as Cox4i1,
# Atp6v1a or Nd-like pseudogenes are not caught by accident.
# -----------------------------------------------------------------------------
MT_CONVENTIONS <- c(
  ensembl   = "^[Mm][Tt]-",
  accession = "^AY172581",
  refseq    = "^(ND[1-6]|ND4L|COX[1-3]|ATP6|ATP8|CYTB)$"
)
MT_PATTERN <- paste(MT_CONVENTIONS, collapse = "|")

add_percent_mt <- function(obj, assay = "RNA", strict = TRUE) {
  genes <- rownames(obj[[assay]])
  hits  <- grep(MT_PATTERN, genes, value = TRUE)
  if (length(hits) == 0) {
    msg <- paste0("No mitochondrial genes matched any known naming convention.\n",
                  "  tried: ",
                  paste(sprintf("%s = %s", names(MT_CONVENTIONS), MT_CONVENTIONS),
                        collapse = " | "), "\n",
                  "  percent.mt would be 0 for every cell and the QC filter would ",
                  "silently do nothing. Inspect rownames() and extend MT_CONVENTIONS.")
    if (strict) stop(msg, call. = FALSE) else warning(msg)
    obj$percent.mt <- 0
    return(obj)
  }
  matched <- names(MT_CONVENTIONS)[vapply(MT_CONVENTIONS,
                                          function(p) any(grepl(p, genes)), logical(1))]
  message(sprintf("percent.mt: %d mito features via '%s' naming (%s%s)",
                  length(hits), paste(matched, collapse = "+"),
                  paste(utils::head(hits, 5), collapse = ", "),
                  if (length(hits) > 5) ", ..." else ""))
  if (length(hits) > 15) {
    warning(length(hits), " features matched the mito pattern -- more than the 13 ",
            "protein-coding mitochondrial genes. Check for false positives before ",
            "trusting percent.mt.")
  }
  PercentageFeatureSet(obj, features = hits, assay = assay, col.name = "percent.mt")
}

# -----------------------------------------------------------------------------
# Misc helpers
# -----------------------------------------------------------------------------
# Metadata columns that must NOT be carried across a CreateSeuratObject() call,
# because they describe the *previous* count matrix.
STALE_META <- function(nms) {
  grepl("^nCount_|^nFeature_|^percent\\.|_snn_res\\.|^seurat_clusters$|^scDblFinder",
        nms)
}

drop_stale_meta <- function(meta) meta[, !STALE_META(names(meta)), drop = FALSE]

# Keep only features actually present, and warn loudly about the rest.
present_features <- function(obj, features, assay = DefaultAssay(obj), label = "features") {
  keep <- intersect(features, rownames(obj[[assay]]))
  missing <- setdiff(features, keep)
  if (length(missing)) {
    message(sprintf("%s: %d/%d not in the %s assay and dropped: %s",
                    label, length(missing), length(features), assay,
                    paste(missing, collapse = ", ")))
  }
  keep
}

# Drop features that reach `min_pct` in NO cluster. dot.min only hides the dot;
# the gene still occupies a column, which is what makes a long neuropeptide
# dotplot unreadable. This removes the column entirely.
expressed_features <- function(obj, features, group.by = "seurat_clusters",
                               min_pct = params$plots$dot_min, assay = "RNA") {
  features <- present_features(obj, features, assay, "features")
  if (!length(features)) return(features)
  cts <- get_expr(obj, assay, "counts")[features, , drop = FALSE]
  grp <- factor(obj@meta.data[[group.by]])
  pct <- vapply(levels(grp), function(g) {
    Matrix::rowMeans(cts[, grp == g, drop = FALSE] > 0)
  }, numeric(length(features)))
  if (is.null(dim(pct))) pct <- matrix(pct, nrow = length(features))
  keep <- features[apply(pct, 1, max) >= min_pct]
  if (length(keep) < length(features)) {
    message(sprintf("expressed_features: %d/%d dropped (below %.0f%% in every cluster): %s",
                    length(features) - length(keep), length(features), 100 * min_pct,
                    paste(setdiff(features, keep), collapse = ", ")))
  }
  keep
}

# -----------------------------------------------------------------------------
# X-axis gene labels for DotPlots
# Seurat's RotatedAxis() is fixed at +45 degrees (labels slope UP to the right).
# A negative angle slopes them DOWN to the right, which reads better under a
# long gene list. The anchoring has to follow the sign or the labels drift away
# from their ticks: hjust = 1 for positive angles, hjust = 0 for negative.
# -----------------------------------------------------------------------------
rotated_x <- function(angle = params$plots$axis_angle,
                      size  = params$plots$axis_text_size) {
  if (is.null(angle)) angle <- 45
  if (is.null(size))  size  <- 8
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = angle,
      hjust = if (angle < 0) 0 else 1,
      vjust = 1,
      size  = size))
}

# UMAP column names differ between Seurat versions (UMAP_1 vs umap_1).
umap_cols <- function(obj, reduction = "umap") colnames(Embeddings(obj, reduction))

# -----------------------------------------------------------------------------
# Memory
# `scale.data` is a DENSE matrix and is by far the largest part of an SCT or
# integrated assay. It is needed only to compute the PCA; DotPlot, VlnPlot and
# FindMarkers all read `data`. Dropping it before saveRDS typically cuts the
# object (and the RAM needed to load it later) by more than half.
# -----------------------------------------------------------------------------
trim_scale_data <- function(obj, assays = c("SCT", "integrated")) {
  before <- as.numeric(utils::object.size(obj))
  for (a in intersect(assays, names(obj@assays))) {
    ok <- tryCatch({
      if (.seurat_v5) {
        if ("scale.data" %in% SeuratObject::Layers(obj[[a]]))
          SeuratObject::LayerData(obj, assay = a, layer = "scale.data") <- NULL
      } else {
        obj <- SetAssayData(obj, assay = a, slot = "scale.data", new.data = new("matrix"))
      }
      TRUE
    }, error = function(e) FALSE)
    if (!ok) message("[mem] could not drop scale.data from assay '", a, "'")
  }
  after <- as.numeric(utils::object.size(obj))
  message(sprintf("[mem] object %.1f -> %.1f GB after dropping scale.data",
                  before / 1024^3, after / 1024^3))
  gc(verbose = FALSE)
  obj
}

log_step <- function(step, obj, note = "") {
  message(sprintf("[%s] %d cells x %d features %s",
                  step, ncol(obj), nrow(obj), note))
  invisible(obj)
}

writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_00.txt")
