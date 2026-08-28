# =============================================================================
# 04_doublets_and_recluster.R  -- doublet detection, removal, and reclustering
#
# Merged from the former 04_doublets_scDbl.R and 04b_recluster_after_doublets.R.
# They were always one logical step: clustering computed WITH doublets present
# is only provisional, and leaving it in place produces hollow clusters (one
# went from 246 cells to 0, another from 128 to 1). Splitting the step across
# two files let the cluster IDs of the two halves drift apart.
#
# What this script does, in order:
#   1. call doublets on RAW RNA counts (not SCT-corrected counts)
#   2. write the pre-removal cluster diagnostics, keyed to the step-03 clusters
#   3. remove doublet CELLS
#   4. recluster the singlets on the existing `integrated` assay
#   5. write the post-reclustering diagnostics, keyed to the FINAL cluster IDs
#
# Step 5 matters: config/clusters_to_drop.csv in step 05 is keyed to these final
# IDs, so this is the table step 05 joins its exclusion reasons against. The
# pre-removal table uses the step-03 numbering and is kept only as a record.
#
# Scope of the reclustering: PCA/UMAP/clusters are recomputed on the existing
# `integrated` assay, NOT re-integrated from scratch. The anchors in step 03
# were computed with doublets present, distributed across samples, which does
# not meaningfully bias the correction. Re-running 03 end to end is the purist
# option and is the right call if the doublet rate is high; state which one you
# used in the Methods.
#
# There is deliberately no checkpoint here. scDblFinder is the slow step, but
# doublets$seed makes it reproducible, so caching its output buys nothing --
# while a cached file that survives a change upstream silently overrides the
# fresh computation, which is exactly how a stale intermediate once propagated
# a wrong clustering through every downstream figure. The script always runs
# from 03_soupx_sct.rds.
# =============================================================================
source("scripts/00_setup.R")

obj <- readRDS("data/processed/03_soupx_sct.rds")
obj <- join_layers_if_split(obj, "RNA")   # Seurat v5 keeps one RNA layer per sample

# --- doublet calling on raw counts -------------------------------------------
DefaultAssay(obj) <- "RNA"
sce <- SingleCellExperiment(list(counts = get_expr(obj, "RNA", "counts")))
colData(sce)$sample <- obj$orig.ident

dims_use <- params$doublets$dims
# scDblFinder simulates artificial doublets, so it is stochastic. Without a
# seed a full re-run marks a slightly different set of cells, which shifts the
# clustering downstream and can dissolve a small cluster that sat near the
# resolution limit. set.seed() alone is not enough -- the function does its own
# sampling, so the seed is passed to it.
if (!is.null(params$doublets$seed)) set.seed(params$doublets$seed)
sce <- scDblFinder(sce, samples = "sample",
                   BPPARAM = BiocParallel::SerialParam(
                     RNGseed = params$doublets$seed),
                   dims = dims_use,
                   includePCs = min(params$doublets$include_pcs, dims_use))

obj$scDblFinder.class <- as.character(sce$scDblFinder.class[match(colnames(obj), colnames(sce))])
obj$scDblFinder.score <- sce$scDblFinder.score[match(colnames(obj), colnames(sce))]
stopifnot(!anyNA(obj$scDblFinder.class))

message(sprintf("[04] global doublet rate: %.1f%%",
                100 * mean(obj$scDblFinder.class == "doublet")))

write_plot(DimPlot(obj, reduction = "umap", group.by = "scDblFinder.class") + plot_theme(),
           "04_umap_doublet-vs-singlet", pdf = TRUE)

n_before        <- ncol(obj)
clusters_before <- table(obj$seurat_clusters)

# =============================================================================
# Cluster-level diagnostics, BEFORE doublet removal
# -----------------------------------------------------------------------------
# This is the evidence a cluster-level exclusion has to rest on. A cluster is
# a plausible doublet artefact when SEVERAL of these agree:
#   - doublet_pct  >> the global rate
#   - median nCount / nFeature far from the median of its own compartment
#   - a large fraction of its cells co-express markers of two lineages that
#     cannot both be true in one nucleus (e.g. Snap25 and Plp1)
#   - no private marker gene: its top genes are a blend of two other clusters
# Being restricted to one sample is NOT evidence of a doublet artefact -- in a
# Cont/Affi/Aggr design that may be the biological result.
#
# These IDs are the step-03 clusters and become stale at the reclustering
# below; the table is kept as a record of how the doublet call behaved.
# =============================================================================
lineage_markers <- list(
  neuron = c("Snap25", "Rbfox3", "Syt1"),
  oligo  = c("Plp1", "Mbp", "Mog"),
  astro  = c("Aqp4", "Slc1a2"),
  micro  = c("Ptprc", "Csf1r", "Apbb1ip")
)
lineage_markers <- lapply(lineage_markers, function(g) present_features(obj, g, "RNA", "lineage markers"))

cts <- get_expr(obj, "RNA", "counts")
lineage_pos <- as.data.frame(lapply(lineage_markers, function(g) {
  if (!length(g)) return(rep(FALSE, ncol(obj)))
  Matrix::colSums(cts[g, , drop = FALSE] > 0) > 0
}))
rownames(lineage_pos) <- colnames(obj)
n_lineages <- rowSums(lineage_pos)

# Depth must be compared WITHIN compartment. Neurons carry ~3x more RNA than
# oligodendrocytes, so a ratio against the global median flags every neuronal
# cluster and carries no information. Each cluster is assigned a compartment
# first, then compared to the median of its own compartment.
cell_tab <- data.frame(
  cluster       = as.character(obj$seurat_clusters),
  doublet       = obj$scDblFinder.class == "doublet",
  nCount        = obj$nCount_RNA,
  nFeature      = obj$nFeature_RNA,
  neuronal      = lineage_pos$neuron,
  multi_lineage = n_lineages >= 2
)

compartment_of <- cell_tab %>%
  group_by(cluster) %>%
  summarise(neuronal_pct = 100 * mean(neuronal), .groups = "drop") %>%
  mutate(compartment = ifelse(neuronal_pct >= 50, "neuronal", "non-neuronal"))

cell_tab <- dplyr::left_join(cell_tab, compartment_of[, c("cluster", "compartment")],
                             by = "cluster")
# as.numeric() strips the dim attribute: tapply() returns a 1-d array, and
# ref_median[compartment] would inherit it, making nCount_vs_compartment a
# matrix column. write_csv() rejects those ("must not contain list or matrix
# columns"), and depth_outlier / flag_score inherit the same shape.
ref_median <- as.numeric(tapply(cell_tab$nCount, cell_tab$compartment, median))
names(ref_median) <- sort(unique(cell_tab$compartment))

diag_tab <- cell_tab %>%
  group_by(cluster) %>%
  summarise(
    n_cells           = n(),
    compartment       = compartment[1],
    doublet_pct       = round(100 * mean(doublet), 1),
    median_nCount     = median(nCount),
    median_nFeature   = median(nFeature),
    multi_lineage_pct = round(100 * mean(multi_lineage), 1),
    .groups = "drop"
  ) %>%
  mutate(
    nCount_vs_compartment = round(median_nCount / as.numeric(ref_median[compartment]), 2),
    doublet_enriched = doublet_pct >= params$doublets$cluster_doublet_pct,
    depth_outlier    = nCount_vs_compartment >= params$doublets$cluster_ncount_ratio |
      nCount_vs_compartment <= 1 / params$doublets$cluster_ncount_ratio,
    lineage_mixed    = multi_lineage_pct >= params$doublets$cluster_multilineage_pct,
    flag_score       = doublet_enriched + depth_outlier + lineage_mixed
  ) %>%
  arrange(desc(flag_score), desc(doublet_pct))

stopifnot("diag_tab has a matrix/list column" =
            all(vapply(diag_tab, function(x) is.null(dim(x)), logical(1))))
readr::write_csv(diag_tab, "data/processed/qc/04_cluster_doublet_diagnostics_preRemoval.csv")
print(as.data.frame(diag_tab))

flagged <- diag_tab$cluster[diag_tab$flag_score >= 2]
if (length(flagged)) {
  message("[04] pre-removal clusters meeting >=2 of 3 doublet criteria: ",
          paste(flagged, collapse = ", "))
} else {
  message("[04] no pre-removal cluster meets >=2 of 3 doublet criteria.")
}

# --- remove doublet CELLS -----------------------------------------------------
rm(cts, lineage_pos); if (exists("sce")) rm(sce); gc(verbose = FALSE)
obj <- subset(obj, subset = scDblFinder.class == "singlet")
Idents(obj) <- "seurat_clusters"

loss <- data.frame(
  cluster = names(clusters_before),
  before  = as.integer(clusters_before),
  after   = as.integer(table(factor(obj$seurat_clusters, levels = names(clusters_before))))
) %>% mutate(pct_lost = round(100 * (before - after) / pmax(before, 1), 1))
readr::write_csv(loss, "data/processed/qc/04_cluster_cell_loss.csv")

message(sprintf("[04] %d -> %d cells after doublet removal", n_before, ncol(obj)))
write_plot(VlnPlot(obj, features = "nCount_RNA", pt.size = 0) + plot_theme(),
           "04_vln_umi_per-cluster_singlets", pdf = TRUE)

# =============================================================================
# Reclustering of the singlets
# =============================================================================
obj$clusters_preDbl <- obj$seurat_clusters
n_clusters_before <- nlevels(droplevels(factor(obj$clusters_preDbl)))

DefaultAssay(obj) <- "integrated"
npcs <- params$integration$npcs
res  <- params$recluster$resolution

# The call sequence below is reproduced EXACTLY as the published analysis ran
# it. Three details are load-bearing and must not be "tidied":
#
#   * the order is PCA -> UMAP -> FindNeighbors -> FindClusters. Reordering it
#     changes how much RNG each step consumes and therefore the UMAP layout.
#   * FindClusters() is called WITHOUT random.seed, i.e. with Seurat's default
#     of 0. Passing 42 gives a different Louvain initialisation and different
#     clusters -- not merely a different picture.
#   * ScaleData() is conditional. In this pipeline 03 trims scale.data, so it
#     always runs in practice; the branch is kept so the behaviour matches.
#
# trim_scale_data() may have dropped scale.data; RunPCA needs it.
has_scale <- tryCatch({
  sd <- get_expr(obj, "integrated", "scale.data")
  !is.null(sd) && nrow(sd) > 0 && ncol(sd) > 0
}, error = function(e) FALSE)
if (!has_scale) {
  message("[04] scale.data missing from the integrated assay -> running ScaleData()")
  obj <- ScaleData(obj, assay = "integrated", verbose = FALSE)
}

obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:npcs, verbose = FALSE)
obj <- FindClusters(obj, resolution = res, verbose = FALSE)
Idents(obj) <- "seurat_clusters"

n_after <- nlevels(obj$seurat_clusters)
message(sprintf("[04b] %d clusters before -> %d after reclustering (resolution %.2f, %d cells)",
                n_before, n_after, res, ncol(obj)))

sizes <- sort(table(obj$seurat_clusters))
tiny <- names(sizes)[sizes < params$annotation$min_cluster_cells]
if (length(tiny)) {
  message("[04b] still below ", params$annotation$min_cluster_cells, " cells: ",
          paste(sprintf("%s (n=%d)", tiny, sizes[tiny]), collapse = ", "),
          " -- if several remain, lower recluster$resolution.")
} else {
  message("[04b] no cluster below ", params$annotation$min_cluster_cells, " cells.")
}

n_after <- nlevels(obj$seurat_clusters)
message(sprintf("[04] %d clusters before -> %d after reclustering (resolution %.2f, %d cells)",
                n_clusters_before, n_after, res, ncol(obj)))

sizes <- sort(table(obj$seurat_clusters))
tiny <- names(sizes)[sizes < params$annotation$min_cluster_cells]
if (length(tiny)) {
  message("[04] below ", params$annotation$min_cluster_cells, " cells: ",
          paste(sprintf("%s (n=%d)", tiny, sizes[tiny]), collapse = ", "),
          " -- if several remain, lower recluster$resolution.")
} else {
  message("[04] no cluster below ", params$annotation$min_cluster_cells, " cells.")
}

# How the old labels map onto the new ones -- shows where the hollow clusters went.
crosstab <- as.data.frame(table(old = obj$clusters_preDbl, new = obj$seurat_clusters)) %>%
  filter(Freq > 0) %>%
  arrange(as.integer(as.character(old)), desc(Freq))
readr::write_csv(crosstab, "data/processed/qc/04b_old_vs_new_clusters.csv")

readr::write_csv(
  data.frame(cluster = names(table(obj$seurat_clusters)),
             n_cells = as.integer(table(obj$seurat_clusters))),
  "data/processed/qc/04b_cluster_sizes.csv")

write_plot(DimPlot(obj, reduction = "umap", label = TRUE) + plot_theme(),
           "04_umap_reclustered", pdf = TRUE)
write_plot(DimPlot(obj, reduction = "umap", group.by = "clusters_preDbl", label = TRUE) +
             ggtitle("Clusters before doublet removal") + plot_theme(),
           "04_umap_preDbl_labels", pdf = TRUE)

# =============================================================================
# Cluster-level diagnostics, AFTER reclustering
# -----------------------------------------------------------------------------
# Keyed to the FINAL cluster IDs, which is what config/clusters_to_drop.csv in
# step 05 refers to. `doublet_pct` is meaningless here (the doublets are gone),
# so the retained cells' mean scDblFinder score is reported instead: a cluster
# built from cells that were only just below the calling threshold still scores
# high. Lineage mixing and sequencing depth carry over unchanged.
# =============================================================================
cts <- get_expr(obj, "RNA", "counts")
lineage_pos <- as.data.frame(lapply(lineage_markers, function(g) {
  if (!length(g)) return(rep(FALSE, ncol(obj)))
  Matrix::colSums(cts[g, , drop = FALSE] > 0) > 0
}))
rownames(lineage_pos) <- colnames(obj)

cell_tab2 <- data.frame(
  cluster       = as.character(obj$seurat_clusters),
  dbl_score     = obj$scDblFinder.score,
  nCount        = obj$nCount_RNA,
  nFeature      = obj$nFeature_RNA,
  neuronal      = lineage_pos$neuron,
  multi_lineage = rowSums(lineage_pos) >= 2
)
compartment2 <- cell_tab2 %>%
  group_by(cluster) %>%
  summarise(neuronal_pct = 100 * mean(neuronal), .groups = "drop") %>%
  mutate(compartment = ifelse(neuronal_pct >= 50, "neuronal", "non-neuronal"))
cell_tab2 <- dplyr::left_join(cell_tab2, compartment2[, c("cluster", "compartment")],
                              by = "cluster")
ref_median2 <- as.numeric(tapply(cell_tab2$nCount, cell_tab2$compartment, median))
names(ref_median2) <- sort(unique(cell_tab2$compartment))

diag_final <- cell_tab2 %>%
  group_by(cluster) %>%
  summarise(
    n_cells             = n(),
    compartment         = compartment[1],
    mean_doublet_score  = round(mean(dbl_score), 3),
    median_nCount       = median(nCount),
    median_nFeature     = median(nFeature),
    multi_lineage_pct   = round(100 * mean(multi_lineage), 1),
    .groups = "drop"
  ) %>%
  mutate(
    nCount_vs_compartment = round(median_nCount / as.numeric(ref_median2[compartment]), 2),
    depth_outlier   = nCount_vs_compartment >= params$doublets$cluster_ncount_ratio |
                      nCount_vs_compartment <= 1 / params$doublets$cluster_ncount_ratio,
    lineage_mixed   = multi_lineage_pct >= params$doublets$cluster_multilineage_pct,
    flag_score      = depth_outlier + lineage_mixed
  ) %>%
  arrange(desc(flag_score), desc(multi_lineage_pct))

stopifnot("diag_final has a matrix/list column" =
            all(vapply(diag_final, function(x) is.null(dim(x)), logical(1))))
readr::write_csv(diag_final, "data/processed/qc/04_cluster_diagnostics.csv")
print(as.data.frame(diag_final))
message("[04] qc/04_cluster_diagnostics.csv is keyed to the FINAL cluster IDs -- ",
        "this is the table config/clusters_to_drop.csv and step 05 refer to.")

rm(cts, lineage_pos); gc(verbose = FALSE)
obj <- trim_scale_data(obj)
saveRDS(obj, "data/processed/04_singlets_soupx_sct.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_04.txt")


