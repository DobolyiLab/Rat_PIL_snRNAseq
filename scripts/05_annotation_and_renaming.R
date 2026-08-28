# =============================================================================
# 05_annotation_and_renaming.R  -- cell-type annotation for all cells
#
# The former 05a_annotation_evidence.R is merged in as section 3: it has to run
# AFTER the reclustering in section 2, because that is what defines the cluster
# IDs the rename table is keyed to.
#
# Order of operations, and why it matters:
#   1. documented cluster exclusion   -- keyed to the step-04 FINAL cluster IDs
#   2. rescale / PCA / cluster / UMAP -- this creates NEW cluster IDs
#   3. lineage evidence + suggested names  -- keyed to the new IDs
#   4. renaming from config/cluster_rename.csv -- also keyed to the new IDs
#
# So the two config files answer to different numbering:
#   config/clusters_to_drop.csv  -> step-04 IDs (qc/04_cluster_diagnostics.csv)
#   config/cluster_rename.csv    -> post-section-2 IDs (qc/05a_*_evidence.csv)
# Regenerate cluster_rename.csv whenever section 2 is re-run.
#
# FIXES vs the original version:
#  * the normalization guard actually works (the old max(data) <= 0 test could
#    never fire, so every RNA dotplot was drawn on raw counts)
#  * the object is saved whether or not a rename table exists, so step 06
#    always has an annotated input
#  * cluster factor levels follow the rename table order, not alphabetical
#  * cluster exclusion is driven by config/clusters_to_drop.csv, which must
#    carry a written reason per cluster, and is cross-checked against the
#    step-04 diagnostics. This replaces undocumented manual removal.
# =============================================================================

source("scripts/00_setup.R")

obj <- readRDS("data/processed/04_singlets_soupx_sct.rds")
DefaultAssay(obj) <- "RNA"
obj <- join_layers_if_split(obj, "RNA")
obj <- ensure_rna_normalized(obj)

# =============================================================================
# 1) Documented cluster exclusion (optional)
# config/clusters_to_drop.csv:  cluster_id,reason
# Nothing is dropped unless the file exists. Every dropped cluster must have a
# reason, and the decision is cross-checked against the step-04 diagnostics.
# =============================================================================
if (file.exists("config/clusters_to_drop.csv")) {
  drop <- readr::read_csv("config/clusters_to_drop.csv", show_col_types = FALSE)
  stopifnot(all(c("cluster_id", "reason") %in% names(drop)))
  drop$cluster_id <- as.character(drop$cluster_id)
  stopifnot("every dropped cluster needs a written reason" = all(nzchar(drop$reason)))
  
  unknown <- setdiff(drop$cluster_id, as.character(unique(obj$seurat_clusters)))
  if (length(unknown)) stop("clusters_to_drop.csv lists non-existent clusters: ",
                            paste(unknown, collapse = ", "),
                            "\nThe list is keyed to a different run of step 04. ",
                            "Re-derive it from qc/04_cluster_diagnostics.csv.",
                            call. = FALSE)
  
  # Sizes are reported so a silently renumbered exclusion list is visible: the
  # reason text names a cell count, and it should match.
  sizes_now <- table(as.character(obj$seurat_clusters))
  message("[05] excluding: ",
          paste(sprintf("%s (n=%d)", drop$cluster_id,
                        as.integer(sizes_now[drop$cluster_id])), collapse = ", "),
          " -- check these counts against the reasons in clusters_to_drop.csv")
  
  # The exclusion list is keyed to the cluster IDs that step 04 produced, so it
  # is joined against qc/04_cluster_diagnostics.csv, which uses the same IDs.
  #
  # It is deliberately NOT joined against qc/05a_cluster_annotation_evidence.csv:
  # that table is written further down in THIS script, after section 2 has
  # reclustered, so its IDs refer to different clusters. Joining it would put
  # another cluster's numbers next to the stated reason, in the very file cited
  # in the Methods.
  diag_path <- "data/processed/qc/04_cluster_diagnostics.csv"
  if (file.exists(diag_path)) {
    dg <- readr::read_csv(diag_path, show_col_types = FALSE)
    dg$cluster <- as.character(dg$cluster)
    unsupported <- drop$cluster_id[drop$cluster_id %in% dg$cluster[dg$flag_score < 1]]
    if (length(unsupported)) {
      warning("These clusters are being dropped although the step-04 diagnostics ",
              "flag nothing about them: ", paste(unsupported, collapse = ", "),
              ". Make sure the stated reason is defensible on its own.")
    }
    drop <- dplyr::left_join(drop, dg, by = c("cluster_id" = "cluster"))
  } else {
    warning(diag_path, " not found -- the exclusion will be recorded without ",
            "its supporting numbers. Re-run 04 to regenerate it.")
  }
  readr::write_csv(drop, "data/processed/qc/05_dropped_clusters_with_evidence.csv")
  
  n0 <- ncol(obj)
  Idents(obj) <- "seurat_clusters"   # explicit: drop table is keyed to cluster IDs
  obj <- subset(obj, idents = setdiff(levels(Idents(obj)), drop$cluster_id))
  message(sprintf("[05] dropped %d cluster(s) (%s): %d -> %d cells",
                  nrow(drop), paste(drop$cluster_id, collapse = ", "), n0, ncol(obj)))
} else {
  message("[05] config/clusters_to_drop.csv not present -> no cluster excluded.")
}

# =============================================================================
# 2) Reclustering and new UMAP after cluster exclusion
# -----------------------------------------------------------------------------
# Removing whole clusters changes which cells define the embedding, so PCA,
# neighbours, clusters and UMAP are recomputed. This creates NEW cluster IDs:
# everything downstream (the evidence table in section 3, cluster_rename.csv,
# every figure) is keyed to these, not to the step-04 numbering.
#
# The parameters below reproduce exactly the values this analysis was run with;
# they live in params.yml so the figures can be regenerated without editing code.
# =============================================================================
obj$cluster_before_reclustering <- as.character(obj$seurat_clusters)

rc <- params$annotation$recluster

# Which assay the reclustering runs on. This is the most consequential setting
# in the script: a different assay gives a different PCA, therefore different
# clusters, therefore different cluster IDs -- and config/cluster_rename.csv is
# keyed to those IDs, so changing it silently relabels every cell type.
#
# "integrated" is what the published analysis used: the batch-corrected assay,
# so the clusters are not driven by sample-of-origin. It falls back to "SCT"
# when no integrated assay is present.
#
# If you ever change this, regenerate config/cluster_rename.csv afterwards --
# the validation further down will refuse to apply the old one.
rc_assay <- if (identical(rc$assay, "auto")) {
  if ("integrated" %in% names(obj@assays)) "integrated" else "SCT"
} else rc$assay
stopifnot("recluster assay not present in the object" = rc_assay %in% names(obj@assays))
DefaultAssay(obj) <- rc_assay
message("[05] reclustering on the '", rc_assay, "' assay ",
        "(config: annotation$recluster$assay)")

set.seed(rc$seed)

# Rescaling is required after subsetting: the previous scaling was centred on
# a cell population that no longer exists.
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = rc$npcs, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:rc$dims)
obj <- FindClusters(obj, resolution = rc$resolution, random.seed = rc$seed)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:rc$dims,
               min.dist = rc$min_dist, seed.use = rc$seed, verbose = FALSE)

message(sprintf("[05] %d clusters after exclusion and reclustering (%d cells)",
                nlevels(obj$seurat_clusters), ncol(obj)))

# Back to RNA for every marker-based operation and figure below.
DefaultAssay(obj) <- "RNA"

markers <- readr::read_csv("config/markers_used_in_figures.csv", show_col_types = FALSE)$gene
markers <- present_features(obj, markers, "RNA", "figure markers")

write_plot(DotPlot(obj, features = markers, group.by = "seurat_clusters",
                   cols = c("lightgrey", "midnightblue"), col.min = -1, scale.min = 0) +
             RotatedAxis() + plot_theme() +
             theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1, size = 14),
                   axis.text.y = element_text(size = 14)),
           "05_dotplot_markers_all", pdf = TRUE)
write_plot(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + plot_theme(),
           "05_umap_by_cluster_all", pdf = TRUE)
write_plot(DimPlot(obj, reduction = "umap", group.by = "orig.ident") + plot_theme(),
           "05_umap_by_sample_all", pdf = TRUE)

# =============================================================================
# 3) Lineage evidence for filling config/cluster_rename.csv
#    (merged in from the former 05a_annotation_evidence.R)
#
# Cluster IDs are NOT stable across pipeline runs. A rename table from an
# earlier run cannot be reused: cluster 3 today is not cluster 3 from before.
# This script scores every cluster against canonical lineage panels, writes the
# evidence to a table, and pre-fills the rename template with a SUGGESTION.
#
# The suggestion is a starting point, not an annotation. Check it against the
# marker table and the dotplot before saving cluster_rename.csv.
# =============================================================================
DefaultAssay(obj) <- "RNA"
obj <- join_layers_if_split(obj, "RNA")
obj <- ensure_rna_normalized(obj)
Idents(obj) <- "seurat_clusters"

# Panels use the same genes as markers_used_in_figures.csv, grouped by lineage.
panels <- list(
  Neuron_pan   = c("Snap25", "Rbfox3", "Syt1"),
  GABA         = c("Gad1", "Gad2", "Slc32a1"),
  GLU          = c("Slc17a6", "Slc17a7"),
  DOPA         = c("Th", "Slc6a3", "Ddc"),
  Astro        = c("Slc1a2", "Aqp4", "Gja1"),
  Oligo        = c("Plp1", "Mbp", "Opalin", "Mog"),
  Oligo_premyel = c("Gpr17", "Bmp4", "Bcas1"),
  OPC          = c("Pdgfra", "Cspg4"),
  Micro        = c("Ptprc", "Apbb1ip", "Csf1r"),
  Endothel     = c("Pecam1", "Flt1", "Cldn5"),
  Pericyte     = c("Abcc9", "Rgs5", "Kcnj8"),
  Ependymal    = c("Foxj1", "Ttr", "Hdc"),
  VLMC         = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Slc6a13", "Ptgds")
)
panels <- lapply(panels, function(g) present_features(obj, g, "RNA", "panel"))
panels <- panels[lengths(panels) > 0]

# --- percent of cells expressing each panel, per cluster ---------------------
cts <- get_expr(obj, "RNA", "counts")
cl  <- factor(obj$seurat_clusters)

pct_expressing <- function(genes) {
  pos <- Matrix::colMeans(cts[genes, , drop = FALSE] > 0)   # fraction of panel genes detected
  tapply(pos, cl, mean) * 100
}
score_mat <- vapply(panels, pct_expressing, numeric(nlevels(cl)))
rownames(score_mat) <- levels(cl)

evidence <- as.data.frame(round(score_mat, 1))
evidence$cluster <- rownames(evidence)
evidence$n_cells <- as.integer(table(cl))

# --- suggested label ---------------------------------------------------------
glia_cols <- intersect(c("Astro","Oligo","Oligo_premyel","OPC","Micro",
                         "Endothel","Pericyte","Ependymal","VLMC"), colnames(score_mat))
neuro_sub <- intersect(c("GABA","GLU","DOPA"), colnames(score_mat))

suggest_one <- function(i) {
  s <- score_mat[i, ]
  is_neuron <- "Neuron_pan" %in% names(s) && s[["Neuron_pan"]] >= 50
  best_glia <- if (length(glia_cols)) glia_cols[which.max(s[glia_cols])] else NA
  glia_val  <- if (length(glia_cols)) max(s[glia_cols]) else 0
  if (is_neuron && glia_val < 40) {
    if (!length(neuro_sub)) return("Neuron")
    return(neuro_sub[which.max(s[neuro_sub])])
  }
  if (glia_val >= 30) return(best_glia)
  "UNRESOLVED"
}
evidence$suggested <- vapply(seq_len(nrow(score_mat)), suggest_one, character(1))

# number clusters within a lineage by size, largest first
evidence <- evidence %>%
  group_by(suggested) %>%
  arrange(desc(n_cells), .by_group = TRUE) %>%
  mutate(suggested_name = if (n() > 1) paste0(suggested, "_", row_number()) else suggested) %>%
  ungroup() %>%
  arrange(as.integer(cluster))

# --- top marker genes per cluster, as a sanity check -------------------------
rm(cts); gc(verbose = FALSE)   # the full count matrix is no longer needed
mk <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 0.5)
readr::write_csv(mk, "data/processed/qc/05a_allcell_markers.csv")

top10 <- mk %>%
  filter(!grepl("^LOC|^AABR|^ENSRNOG", gene)) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop") %>%
  mutate(cluster = as.character(cluster))

evidence <- dplyr::left_join(evidence, top10, by = "cluster") %>%
  select(cluster, n_cells, suggested_name, everything())

print(as.data.frame(evidence[, c("cluster", "n_cells", "suggested_name", "Neuron_pan",
                                 intersect(c("GABA","GLU","DOPA","Astro","Oligo","OPC","Micro"),
                                           colnames(evidence)))]))

# --- flags: mixed lineage and tiny clusters ----------------------------------
# A cluster that scores high on BOTH a neuronal and a glial panel, or on two
# incompatible glial panels, is a doublet/ambient candidate -- the label from
# suggest_one() alone hides this.
excl <- setdiff(colnames(score_mat), "Neuron_pan")
flags <- t(apply(score_mat, 1, function(s) {
  o <- sort(s[excl], decreasing = TRUE)
  c(top_lineage = names(o)[1], top_score = unname(o[1]),
    second_lineage = names(o)[2], second_score = unname(o[2]))
}))
flags <- as.data.frame(flags, stringsAsFactors = FALSE)
flags$cluster <- rownames(score_mat)
flags$top_score    <- as.numeric(flags$top_score)
flags$second_score <- as.numeric(flags$second_score)

neuronal_panels <- intersect(c("GABA", "GLU", "DOPA"), excl)
is_neuronal <- function(x) x %in% neuronal_panels

# Panels that describe ONE differentiation continuum or ONE compartment are not
# mutually exclusive: an OPC legitimately scores on Oligo_premyel (Gpr17), and a
# mural cell scores on both Pericyte and Endothel. Co-scoring within a family is
# biology, not a doublet, so those pairs must not raise flag_mixed.
FAMILIES <- list(oligo_lineage = c("OPC", "Oligo_premyel", "Oligo"),
                 vascular      = c("Endothel", "Pericyte", "VLMC"))
same_family <- function(a, b) {
  any(vapply(FAMILIES, function(f) a %in% f && b %in% f, logical(1)))
}

flags$flag_mixed <- mapply(function(a, b, sa, sb) {
  if (same_family(a, b)) return(FALSE)
  sb >= params$annotation$mixed_second_min &&
    sb >= params$annotation$mixed_ratio * sa &&
    (is_neuronal(a) != is_neuronal(b) || (!is_neuronal(a) && !is_neuronal(b)))
}, flags$top_lineage, flags$second_lineage, flags$top_score, flags$second_score)

evidence <- dplyr::left_join(
  evidence,
  flags[, c("cluster", "top_lineage", "top_score", "second_lineage",
            "second_score", "flag_mixed")],
  by = "cluster")
evidence$flag_tiny <- evidence$n_cells < params$annotation$min_cluster_cells

if (any(evidence$flag_mixed)) {
  message("[05a] MIXED lineage (doublet/ambient candidates): ",
          paste(sprintf("%s (%s %.0f%% + %s %.0f%%)",
                        evidence$cluster[evidence$flag_mixed],
                        evidence$top_lineage[evidence$flag_mixed],
                        evidence$top_score[evidence$flag_mixed],
                        evidence$second_lineage[evidence$flag_mixed],
                        evidence$second_score[evidence$flag_mixed]),
                collapse = "; "))
}
if (any(evidence$flag_tiny)) {
  message("[05a] clusters below ", params$annotation$min_cluster_cells, " cells: ",
          paste(sprintf("%s (n=%d)", evidence$cluster[evidence$flag_tiny],
                        evidence$n_cells[evidence$flag_tiny]), collapse = ", "),
          " -- consider lowering integration$resolution rather than annotating these.")
}

readr::write_csv(evidence, "data/processed/qc/05a_cluster_annotation_evidence.csv")

unresolved <- evidence$cluster[grepl("^UNRESOLVED", evidence$suggested_name)]
if (length(unresolved)) {
  message("[05a] no confident lineage for cluster(s): ", paste(unresolved, collapse = ", "),
          " -- inspect their top markers and the step-04 doublet diagnostics.")
}

# --- pre-filled template -----------------------------------------------------
# Built directly, without rename(), which can resolve to S4Vectors::rename().
# n_cells is written alongside the names: it is the fingerprint the renaming
# step below uses to detect a rename table left over from a different run.
suggested_tbl <- data.frame(cluster_id = evidence$cluster,
                            new_name   = evidence$suggested_name,
                            n_cells    = evidence$n_cells,
                            stringsAsFactors = FALSE)
readr::write_csv(suggested_tbl, "config/cluster_rename_suggested.csv")
message("[05a] wrote config/cluster_rename_suggested.csv -- REVIEW IT, edit the names, ",
        "then save it as config/cluster_rename.csv and re-run 05.")

# --- dotplot grouped by lineage ---------------------------------------------
write_plot(DotPlot(obj, features = panels, group.by = "seurat_clusters",
                   cols = c("lightgrey", "midnightblue"), dot.scale = 6) +
             RotatedAxis() + plot_theme() +
             theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1, size = 8),
                   axis.text.y = element_text(size = 14)),
           "05a_dotplot_lineage_panels", pdf = TRUE)

# =============================================================================
# 4) Renaming
# =============================================================================
cur <- as.character(sort(unique(as.integer(as.character(obj$seurat_clusters)))))

# -----------------------------------------------------------------------------
# Which table names the clusters.
#
#   config/cluster_rename.csv            a table YOU wrote. Takes precedence.
#   config/cluster_rename_suggested.csv  written by section 3 above, from the
#                                        lineage panels, against the clusters
#                                        that exist right now.
#
# The suggested table is regenerated on every run, so it can never be stale.
# A hand-written cluster_rename.csv can be, which is why the validation below
# refuses to apply one that does not match the current clustering.
#
# No template is written: section 3 already produces a filled-in starting point.
# -----------------------------------------------------------------------------
rename_path <- if (file.exists("config/cluster_rename.csv")) {
  "config/cluster_rename.csv"
} else if (file.exists("config/cluster_rename_suggested.csv")) {
  "config/cluster_rename_suggested.csv"
} else {
  NA_character_
}

if (is.na(rename_path)) {
  
  message("[05] no rename table found -- clusters keep their numeric labels. ",
          "Section 3 normally writes config/cluster_rename_suggested.csv; if it ",
          "is missing, check that section for errors.")
  obj$celltype_detailed <- factor(obj$seurat_clusters, levels = cur)
  
} else {
  
  message("[05] naming clusters from ", rename_path,
          if (basename(rename_path) == "cluster_rename_suggested.csv")
            " (auto-generated; write config/cluster_rename.csv to override)" else
              " (hand-written; it overrides the suggestion)")
  
  ren <- readr::read_delim(
    rename_path,
    delim = NULL,
    show_col_types = FALSE
  )
  stopifnot(all(c("cluster_id", "new_name") %in% names(ren)))
  ren$cluster_id <- as.character(ren$cluster_id)
  ren <- ren[ren$cluster_id %in% cur & nzchar(ren$new_name), , drop = FALSE]
  
  # -------------------------------------------------------------------------
  # Validation. Cluster IDs are not stable across runs: any change upstream
  # (step 04, the exclusion list, the reclustering parameters above) renumbers
  # them, and a rename table from an earlier run then puts the right names on
  # the wrong clusters. That failure is silent on every table and only visible
  # as a biologically impossible UMAP -- GABAergic neurons sharing a blob with
  # oligodendrocytes, and a legend with bare numbers in it.
  #
  # So: refuse to apply a table that does not cover the current clusters.
  # -------------------------------------------------------------------------
  unmapped <- setdiff(cur, ren$cluster_id)
  extra    <- setdiff(ren$cluster_id, cur)
  
  size_mismatch <- character(0)
  if ("n_cells" %in% names(ren)) {
    now <- table(as.character(obj$seurat_clusters))
    common <- intersect(ren$cluster_id, names(now))
    expected <- setNames(ren$n_cells[match(common, ren$cluster_id)], common)
    size_mismatch <- common[as.integer(now[common]) != as.integer(expected)]
  }
  
  if (length(unmapped) || length(extra) || length(size_mismatch)) {
    msg <- c("config/cluster_rename.csv does not match the current clustering.",
             if (length(unmapped))
               paste0("  clusters with no name: ", paste(unmapped, collapse = ", ")),
             if (length(extra))
               paste0("  names for clusters that no longer exist: ",
                      paste(extra, collapse = ", ")),
             if (length(size_mismatch))
               paste0("  cluster sizes changed since the table was written: ",
                      paste(size_mismatch, collapse = ", ")),
             "",
             "The table was built for a different run. Applying it would label",
             "cells with another cluster's cell type.",
             "",
             "To recover:",
             "  1. review config/cluster_rename_suggested.csv, just written by",
             "     section 3 against the CURRENT clusters",
             "  2. edit the names, save it as config/cluster_rename.csv",
             "  3. re-run this script",
             "",
             "Set annotation$allow_partial_rename: true in params.yml to proceed",
             "anyway, leaving unmatched clusters with numeric labels.")
    if (isTRUE(params$annotation$allow_partial_rename)) {
      warning(paste(msg, collapse = "\n"))
    } else {
      stop(paste(msg, collapse = "\n"), call. = FALSE)
    }
  }
  
  map <- setNames(ren$new_name, ren$cluster_id)
  lab <- ifelse(as.character(obj$seurat_clusters) %in% names(map),
                map[as.character(obj$seurat_clusters)],
                as.character(obj$seurat_clusters))
  # level order follows the rename table, then any unmapped numeric clusters
  obj$celltype_detailed <- factor(lab, levels = c(unique(ren$new_name), unmapped))
  Idents(obj) <- "celltype_detailed"
  
  # (05_umap_celltype_labeled_all is written below, in its styled form)
  write_plot(DotPlot(obj, features = markers, group.by = "celltype_detailed",
                     cols = c("lightgrey", "midnightblue"), col.min = -1, scale.min = 0) +
               RotatedAxis() + plot_theme() +
               theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1, size = 14)),
             "05_dotplot_by_celltype_all", pdf = TRUE)
}

write_plot(DimPlot(obj, reduction = "umap", label = TRUE,
                   label.size = 4.5, repel = TRUE) + plot_theme() +
             theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1, size = 14)),
           "05_umap_celltype_labeled_all", pdf = TRUE)

p <- DimPlot(
  obj,
  reduction = "umap",
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  plot_theme() +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 3)
    )
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

ggsave(
  "figures/05_umap_celltype_labeled_all.tiff",
  plot = p,
  width = 12,
  height = 8,
  units = "in",
  dpi = 500,
  compression = "lzw"
)

readr::write_csv(
  as.data.frame(table(celltype = obj$celltype_detailed, sample = obj$orig.ident)),
  "data/processed/qc/05_celltype_by_sample_counts.csv")

saveRDS(obj, "data/processed/05_annotated_all.rds")   # always saved
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_05.txt")