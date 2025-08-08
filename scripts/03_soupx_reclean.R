source("scripts/00_setup.R")
sct <- readRDS("data/processed/02_sct_integrated.rds")
sct$barcode_short <- sub("_\\d+$", "", colnames(sct))
get_cluster_map <- function(obj, sample_id) {
  meta <- obj@meta.data; meta <- meta[meta$orig.ident == sample_id, , drop = FALSE]
  bar_short <- meta$barcode_short; names(bar_short) <- rownames(meta)
  cl <- meta$seurat_clusters; names(cl) <- bar_short
  stats::setNames(cl, bar_short)
}
adjusted_list <- vector("list", length(params$samples))
seu_list <- vector("list", length(params$samples))
for (i in seq_along(params$samples)) {
  sample_id <- params$samples[[i]]$id; outs_path <- params$samples[[i]]$path
  filt <- Read10X(file.path(outs_path, "filtered_feature_bc_matrix"))
  raw  <- Read10X(file.path(outs_path, "raw_feature_bc_matrix"))

  sc  <- SoupChannel(tod = raw, toc = filt)
  sc$metaData$cell <- rownames(sc$metaData)
  cluster_map <- get_cluster_map(sct, sample_id)
  sc$metaData$cluster <- as.character(cluster_map[sc$metaData$cell])
  sc <- setClusters(sc, clusters = sc$metaData$cluster)
  sc <- autoEstCont(sc, forceAccept = TRUE)
  out <- adjustCounts(sc, roundToInt = TRUE)
  adjusted_list[[i]] <- out
  keep_cells <- intersect(colnames(out), sct$barcode_short[sct$orig.ident == sample_id])
  meta_sub <- sct@meta.data[sct$orig.ident == sample_id, , drop = FALSE]
  rownames(meta_sub) <- meta_sub$barcode_short
  meta_keep <- meta_sub[keep_cells, , drop = FALSE]
  so <- CreateSeuratObject(out[, keep_cells, drop = FALSE], meta.data = meta_keep)
  seu_list[[i]] <- so
}
soupx <- seu_list
soupx <- lapply(soupx, function(x){ PercentageFeatureSet(x, pattern="^mt-|^MT-", col.name="percent.mt") })
soupx <- lapply(soupx, function(x){ subset(x, subset = nFeature_RNA > params$qc$min_features & nFeature_RNA < params$qc$max_features & percent.mt < params$qc$max_percent_mt) })
soupx <- lapply(soupx, function(x){ SCTransform(x, vst.flavor="v2", vars.to.regress="nCount_RNA", verbose=FALSE) %>% RunPCA(npcs = params$integration$npcs, verbose=FALSE) })
features <- SelectIntegrationFeatures(object.list = soupx, nfeatures = 3000)
soupx <- PrepSCTIntegration(object.list = soupx, anchor.features = features)
soupx.anchors <- FindIntegrationAnchors(object.list = soupx, normalization.method = "SCT", anchor.features = features)
soupx.sct <- IntegrateData(anchorset = soupx.anchors, normalization.method = "SCT")
soupx.sct <- RunPCA(soupx.sct, npcs = params$integration$npcs, verbose = FALSE)
soupx.sct <- RunUMAP(soupx.sct, reduction = "pca", dims = 1:params$integration$npcs, verbose = FALSE)
soupx.sct <- FindNeighbors(soupx.sct, reduction = "pca", dims = 1:params$integration$npcs)
soupx.sct <- FindClusters(soupx.sct, resolution = params$integration$resolution)
p_umap <- DimPlot(soupx.sct, reduction="umap", group.by="ident", label=TRUE) + plot_theme()
write_plot(p_umap, "03_umap_soupx_reintegrated", pdf=TRUE)
p_vln <- VlnPlot(soupx.sct, features="nCount_RNA", pt.size=0) + plot_theme()
write_plot(p_vln, "03_vln_nCount_per-cluster", pdf=TRUE)
saveRDS(soupx.sct, "data/processed/03_soupx_sct.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_03.txt")
