# =============================================================================
# 02_integration_SCT.R  -- preliminary SCT integration + clustering
# These clusters exist only to give SoupX its cluster structure in step 03.
# They are stored as `clusters_prelim` so a later FindClusters() cannot
# silently overwrite them.
#
# FIXES vs original:
#  * RNA assay normalized before any RNA-based DotPlot (was raw counts)
#  * vars.to.regress = "nCount_RNA" removed: SCTransform v2 already models
#    sequencing depth; regressing it again can remove real biological variance
#  * nfeatures / dims taken from params.yml so 02 and 03 cannot drift apart
#  * FindIntegrationAnchors dims tied to npcs
#  * clusters saved under an explicit name
# =============================================================================
source("scripts/00_setup.R")

combined <- readRDS("data/processed/01_combined_qc.rds")
npcs  <- params$integration$npcs
nfeat <- params$integration$nfeatures

lst <- SplitObject(combined, split.by = "orig.ident")
lst <- lapply(lst, function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
  RunPCA(x, npcs = npcs, verbose = FALSE)
})

features <- SelectIntegrationFeatures(lst, nfeatures = nfeat)
lst      <- PrepSCTIntegration(lst, anchor.features = features)
anchors  <- FindIntegrationAnchors(object.list = lst, normalization.method = "SCT",
                                   anchor.features = features, dims = 1:npcs)
sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:npcs)

DefaultAssay(sct) <- "integrated"
sct <- RunPCA(sct, npcs = npcs, verbose = FALSE)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:npcs, verbose = FALSE)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:npcs, verbose = FALSE)
sct <- FindClusters(sct, resolution = params$integration$resolution, verbose = FALSE)

# Freeze the preliminary labels under their own name.
sct$clusters_prelim <- sct$seurat_clusters
log_step("02", sct, sprintf("(%d preliminary clusters)", nlevels(sct$clusters_prelim)))

write_plot(DimPlot(sct, group.by = "orig.ident", reduction = "umap") + plot_theme(),
           "02_umap_by-sample", pdf = TRUE)
write_plot(DimPlot(sct, group.by = "clusters_prelim", reduction = "umap", label = TRUE) + plot_theme(),
           "02_umap_by-cluster", pdf = TRUE)

# ---- marker dotplot on properly normalized RNA ------------------------------
DefaultAssay(sct) <- "RNA"
sct <- ensure_rna_normalized(sct)

markers <- readr::read_csv("config/markers_used_in_figures.csv", show_col_types = FALSE)$gene
markers <- present_features(sct, markers, assay = "RNA", label = "figure markers")

p_dot <- DotPlot(sct, features = markers, group.by = "clusters_prelim",
                 cols = c("lightgrey", "midnightblue"), col.min = -1, scale.min = 0) +
  RotatedAxis() + plot_theme()
write_plot(p_dot, "02_dotplot_markers_in_paper", pdf = TRUE)

saveRDS(sct, "data/processed/02_sct_integrated.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_02.txt")
