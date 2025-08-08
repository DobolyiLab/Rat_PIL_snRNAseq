source("scripts/00_setup.R")
objs <- lapply(params$samples, function(s) {
  so <- CreateSeuratObject(counts = read_10x(s$path), project = s$id, min.cells = 1, min.features = 1)
  so <- PercentageFeatureSet(so, pattern = "^mt-", col.name = "percent.mt")
  so <- PercentageFeatureSet(so, pattern = "^AY[0-9]", col.name = "percent.ay")
  so$percent.mt <- as.numeric(so$percent.mt) + as.numeric(so$percent.ay)
  subset(so, subset = nFeature_RNA > params$qc$min_features & percent.mt < params$qc$max_percent_mt)
})
combined <- Reduce(function(a,b) merge(a, y=b), objs)
p_vln <- VlnPlot(combined, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                 group.by = "orig.ident", pt.size = 0, ncol = 3) + plot_theme()
write_plot(p_vln, "01_violin_feature-count-mt", pdf = TRUE)
saveRDS(combined, "data/processed/01_combined_qc.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_01.txt")
