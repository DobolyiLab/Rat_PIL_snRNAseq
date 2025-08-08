source("scripts/00_setup.R")
combined <- readRDS("data/processed/01_combined_qc.rds")
lst <- SplitObject(combined, split.by = "orig.ident")
lst <- lapply(lst, function(x) SCTransform(x, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE) %>%
  RunPCA(npcs = params$integration$npcs, verbose = FALSE))
features <- SelectIntegrationFeatures(lst, nfeatures = 2000)
lst <- PrepSCTIntegration(lst, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = lst, normalization.method = "SCT", anchor.features = features)
sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
sct <- RunPCA(sct, npcs = params$integration$npcs)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:params$integration$npcs)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:params$integration$npcs)
sct <- FindClusters(sct, resolution = params$integration$resolution)
p1 <- DimPlot(sct, group.by = "orig.ident", reduction = "umap", label = TRUE) + plot_theme()
p2 <- DimPlot(sct, group.by = "seurat_clusters", reduction = "umap", label = TRUE) + plot_theme()
write_plot(p1, "02_umap_by-sample", pdf = TRUE)
write_plot(p2, "02_umap_by-cluster", pdf = TRUE)
DefaultAssay(sct) <- "RNA"
markers <- readr::read_csv("config/markers_used_in_figures.csv", show_col_types = FALSE)$gene
p_dot <- DotPlot(sct, features = markers, cols = c("lightgrey","midnightblue"),
                 col.min = -1, scale.min = 0) + plot_theme()
write_plot(p_dot, "02_dotplot_markers_in_paper", pdf = TRUE)
saveRDS(sct, "data/processed/02_sct_integrated.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_02.txt")
