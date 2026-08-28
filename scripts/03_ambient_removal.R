# =============================================================================
# 03_ambient_removal.R  -- ambient RNA removal + re-integration
#
# SoupX requires the RAW (unfiltered) droplet matrix to estimate the soup
# profile. Only filtered matrices are available in this dataset, so the method
# is selected automatically:
#
#   soupx   -- raw matrices present (params$samples[[i]]$raw_dir set and valid)
#   decontx -- filtered only; celda::decontX estimates contamination from the
#              expression profiles of the clusters themselves
#   none    -- no ambient correction; the object passes through unchanged
#
# Whichever runs is logged to qc/03_ambient_method.csv and must be stated in
# the Methods. Do not silently switch between them between runs.
#
# Other fixes vs the original: safe barcode mapping, honest rho estimation,
# stale metadata dropped, mito pattern shared with step 01.
# =============================================================================
source("scripts/00_setup.R")

sct  <- readRDS("data/processed/02_sct_integrated.rds")
npcs <- params$integration$npcs
stopifnot("clusters_prelim" %in% colnames(sct@meta.data))

method <- params$ambient$method
if (identical(method, "auto")) {
  method <- if (raw_matrices_available()) "soupx" else params$ambient$method_without_raw
}
message("[03] ambient correction method: ", method)
if (method != "soupx") {
  message("[03] NOTE: SoupX was not used because the raw droplet matrices are ",
          "not available. This must appear in the Methods.")
}

seu_list <- vector("list", length(params$samples))
est_log  <- list()

for (i in seq_along(params$samples)) {
  s <- params$samples[[i]]
  sample_id <- s$id
  
  # --- cluster labels for this sample, keyed by plain barcode ---------------
  meta_sub <- sct@meta.data[sct$orig.ident == sample_id, , drop = FALSE]
  bc <- sub(paste0("^", sample_id, "_"), "", rownames(meta_sub))
  stopifnot(!anyDuplicated(bc))
  rownames(meta_sub) <- bc
  cluster_map <- setNames(as.character(meta_sub$clusters_prelim), bc)
  
  filt <- read_sample_matrix(s, "filtered")
  cells_with_cluster <- intersect(colnames(filt), names(cluster_map))
  message(sprintf("[03/%s] %d/%d filtered cells carry a cluster label",
                  sample_id, length(cells_with_cluster), ncol(filt)))
  stopifnot(length(cells_with_cluster) > 0)
  toc  <- filt[, cells_with_cluster, drop = FALSE]
  clus <- cluster_map[colnames(toc)]
  
  # --- ambient correction ---------------------------------------------------
  if (method == "soupx") {
    raw <- read_sample_matrix(s, "raw")
    
    # The raw and filtered matrices must come from the SAME cellranger run.
    # Two candidate R3 runs exist on disk, so this is checked, not assumed.
    if (!identical(rownames(raw), rownames(filt))) {
      stop(sprintf(paste0("[%s] raw and filtered matrices have different feature sets ",
                          "(raw: %d, filtered: %d). They are not from the same cellranger run -- ",
                          "fix raw_dir in params.yml (see check_raw_match.R)."),
                   sample_id, nrow(raw), ncol(filt)), call. = FALSE)
    }
    if (!all(colnames(toc) %in% colnames(raw))) {
      stop(sprintf(paste0("[%s] %d filtered barcodes are absent from the raw matrix. ",
                          "Wrong raw run selected in params.yml."),
                   sample_id, sum(!colnames(toc) %in% colnames(raw))), call. = FALSE)
    }
    
    sc  <- SoupChannel(tod = raw, toc = toc)
    sc  <- setClusters(sc, clusters = clus)
    stopifnot(!anyNA(sc$metaData$clusters))
    
    forced <- FALSE
    sc <- tryCatch(
      autoEstCont(sc, forceAccept = FALSE),
      error = function(e) {
        warning(sprintf("[03/%s] autoEstCont rejected the estimate (%s); retrying forced.",
                        sample_id, conditionMessage(e)))
        forced <<- TRUE
        autoEstCont(sc, forceAccept = TRUE)
      })
    rho <- unique(sc$metaData$rho)[1]
    out <- adjustCounts(sc, roundToInt = TRUE)
    est_log[[sample_id]] <- data.frame(sample = sample_id, method = "soupx",
                                       contamination = rho, forced = forced)
    message(sprintf("[03/%s] SoupX rho = %.3f%s", sample_id, rho,
                    if (forced) " (forced)" else ""))
    
  } else if (method == "decontx") {
    if (!requireNamespace("celda", quietly = TRUE))
      stop("celda is required for the decontX fallback: BiocManager::install('celda')",
           call. = FALSE)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = toc))
    sce <- celda::decontX(sce, z = as.factor(clus), seed = 42, verbose = FALSE)
    out <- round(celda::decontXcounts(sce))
    dimnames(out) <- dimnames(toc)
    cont <- mean(sce$decontX_contamination)
    est_log[[sample_id]] <- data.frame(sample = sample_id, method = "decontx",
                                       contamination = cont, forced = NA)
    message(sprintf("[03/%s] decontX mean contamination = %.3f", sample_id, cont))
    
  } else {
    out <- toc
    est_log[[sample_id]] <- data.frame(sample = sample_id, method = "none",
                                       contamination = 0, forced = NA)
  }
  
  # --- rebuild the Seurat object -------------------------------------------
  keep_cells  <- colnames(out)
  counts_keep <- out[, keep_cells, drop = FALSE]
  colnames(counts_keep) <- paste0(sample_id, "_", keep_cells)
  
  meta_keep <- drop_stale_meta(meta_sub[keep_cells, , drop = FALSE])
  meta_keep$clusters_prelim <- cluster_map[keep_cells]
  rownames(meta_keep) <- colnames(counts_keep)
  
  seu_list[[i]] <- CreateSeuratObject(counts_keep, project = sample_id, meta.data = meta_keep)
}
names(seu_list) <- SAMPLE_IDS
readr::write_csv(dplyr::bind_rows(est_log), "data/processed/qc/03_ambient_method.csv")

# --- QC again on the corrected counts ----------------------------------------
seu_list <- lapply(seu_list, function(x) {
  x <- add_percent_mt(x, strict = TRUE)
  subset(x, subset = nFeature_RNA > params$qc$min_features &
           nFeature_RNA < params$qc$max_features &
           percent.mt   < params$qc$max_percent_mt)
})
seu_list <- lapply(seu_list, function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
  RunPCA(x, npcs = npcs, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = params$integration$nfeatures)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = features)
anchors  <- FindIntegrationAnchors(object.list = seu_list, normalization.method = "SCT",
                                   anchor.features = features, dims = 1:npcs)
clean.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:npcs)

DefaultAssay(clean.sct) <- "integrated"
clean.sct <- RunPCA(clean.sct, npcs = npcs, verbose = FALSE)
clean.sct <- RunUMAP(clean.sct, reduction = "pca", dims = 1:npcs, verbose = FALSE)
clean.sct <- FindNeighbors(clean.sct, reduction = "pca", dims = 1:npcs, verbose = FALSE)
clean.sct <- FindClusters(clean.sct, resolution = params$integration$resolution, verbose = FALSE)
clean.sct$ambient_method <- method
log_step("03", clean.sct, sprintf("(%d clusters, ambient = %s)",
                                  nlevels(clean.sct$seurat_clusters), method))

write_plot(DimPlot(clean.sct, reduction = "umap", label = TRUE) + plot_theme(),
           "03_umap_reintegrated", pdf = TRUE)
write_plot(VlnPlot(clean.sct, features = "nCount_RNA", pt.size = 0) + plot_theme(),
           "03_vln_nCount_per-cluster", pdf = TRUE)

saveRDS(clean.sct, "data/processed/03_soupx_sct.rds")   # filename kept so 04 keeps working
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_03.txt")
