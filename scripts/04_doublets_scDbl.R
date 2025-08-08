source("scripts/00_setup.R")
obj <- readRDS("data/processed/03_soupx_sct.rds")
DefaultAssay(obj) <- "SCT"
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce, samples = "orig.ident", includePCs = 20)
obj$scDblFinder.class <- sce$scDblFinder.class[match(colnames(obj), colnames(sce))]
Idents(obj) <- "scDblFinder.class"
p <- DimPlot(obj, reduction = "umap", label = TRUE) + plot_theme()
write_plot(p, "04_umap_doublet-vs-singlet", pdf = TRUE)
obj <- subset(obj, idents = "singlet"); Idents(obj) <- "seurat_clusters"
p2 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0) + plot_theme()
write_plot(p2, "04_vln_umi_per-cluster_singlets", pdf = TRUE)
saveRDS(obj, "data/processed/04_singlets_soupx_sct.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_04.txt")
