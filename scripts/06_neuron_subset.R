source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# 06_subneurons_analysis.R
# - subset neurons
# - per-sample proportions
# - subset clusters with >=10% contribution (data-driven)
# - re-SCT/UMAP
# - Fos/Calb1 viz + thresholded violin + UMAP highlighting
# - neurotransmitter + neuropeptide dotplots
# - top markers per cluster (3 genes) dotplot

# 1) Subset neurons
obj <- readRDS("data/processed/04_singlets_soupx_sct.rds")
DefaultAssay(obj) <- "RNA"
expr <- FetchData(obj, vars = c("Snap25","Rbfox3"))
keep <- rownames(expr)[rowSums(expr > 0) > 0]
neurons <- subset(obj, cells = keep)
neurons <- SCTransform(neurons, vst.flavor="v2", vars.to.regress="nCount_RNA", verbose=FALSE)
neurons <- RunPCA(neurons, npcs=50, verbose=FALSE)
neurons <- RunUMAP(neurons, dims=1:20, min.dist=0.5)
neurons <- FindNeighbors(neurons, dims=1:20)
neurons <- FindClusters(neurons, resolution=0.6)
p_umap <- DimPlot(neurons, label=TRUE) + ggtitle("Neuron subclustering") + plot_theme()
write_plot(p_umap, "06_umap_neuron_subclustering", pdf=TRUE)


# 2) Per-sample proportions per cluster
tab <- as.data.frame(table(neurons$seurat_clusters, neurons$orig.ident))
colnames(tab) <- c("cluster","sample","cell_count")
cluster_total <- tab %>% group_by(cluster) %>% summarise(total = sum(cell_count), .groups="drop")
cluster_sample_pct <- tab %>% left_join(cluster_total, by="cluster") %>%
  mutate(percent = 100 * cell_count / total) %>% select(cluster, sample, percent)
readr::write_csv(cluster_sample_pct, "data/processed/06_cluster_sample_percent.csv")

# 3) Choose clusters with >=10% contribution in at least one sample (data-driven)
keep_clusters <- cluster_sample_pct %>% group_by(cluster) %>% summarise(max_pct = max(percent)) %>%
  filter(max_pct >= 10) %>% pull(cluster) %>% as.character()
message("Keeping clusters (>=10% in at least one sample): ", paste(keep_clusters, collapse=", "))

subneurons <- subset(neurons, idents = keep_clusters)

# Re-normalize / recluster
subneurons <- SCTransform(subneurons, vst.flavor="v2", vars.to.regress="nCount_RNA", verbose=FALSE)
subneurons <- RunPCA(subneurons, npcs=50, verbose=FALSE)
subneurons <- RunUMAP(subneurons, dims=1:20)
subneurons <- FindNeighbors(subneurons, dims=1:20)
subneurons <- FindClusters(subneurons, resolution=0.6)

markers <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.5,
                          latent.vars=c("nCount_RNA"))
readr::write_csv(markers, "data/processed/markers_neurons.csv")

# Save UMAP
p_sub <- DimPlot(subneurons, label=TRUE) + ggtitle("Subneurons (>=10% cluster inclusion)") + plot_theme()
write_plot(p_sub, "06_umap_subneurons", pdf=TRUE)

# 4) Fos/Calb1 analysis: thresholding by cluster with >=8 cells > 0.25
DefaultAssay(subneurons) <- "RNA"
genes_fc <- c("Fos","Calb1")
expr_data <- FetchData(subneurons, vars = c(genes_fc, "seurat_clusters"))
expr_data$cell <- rownames(expr_data)

for (g in genes_fc) {
  expr_data <- expr_data %>% group_by(seurat_clusters) %>%
    mutate(!!g := ifelse(sum(.data[[g]] > 0.25) < 8, 0, .data[[g]])) %>%
    ungroup()
}

mat <- t(as.matrix(expr_data[, genes_fc, drop=FALSE]))
colnames(mat) <- expr_data$cell
subneurons[["temp_expr"]] <- CreateAssayObject(counts = mat)

p_vln <- VlnPlot(subneurons, features = genes_fc, assay = "temp_expr", pt.size = 1, flip = TRUE, group.by = "seurat_clusters")
write_plot(p_vln, "06_violin_Fos_Calb1_thresholded", pdf=TRUE)

# UMAP highlight
df_umap <- FetchData(subneurons, vars = c("umap_1","umap_2"))
df_umap$cell <- rownames(df_umap)

marker_status <- expr_data %>%
  pivot_longer(cols = all_of(genes_fc), names_to = "gene", values_to = "expr") %>%
  group_by(seurat_clusters, gene) %>%
  mutate(valid = sum(expr > 0.25) >= 8) %>%
  ungroup() %>%
  mutate(marker_group = ifelse(expr > 0.25 & valid, gene, "None")) %>%
  select(cell, gene, marker_group)

marker_wide <- marker_status %>%
  pivot_wider(names_from = gene, values_from = marker_group) %>%
  mutate(marker_group = case_when(
    Fos != "None" ~ "Fos",
    Calb1 != "None" ~ "Calb1",
    TRUE ~ "None"
  ))

df_umap <- dplyr::left_join(df_umap, marker_wide[, c("cell","marker_group")], by="cell")

tiff("figures/06_umap_Fos_Calb1.tiff", width = 3000, height = 2000, res = 300)
print(ggplot(df_umap, aes(x = umap_1, y = umap_2)) +
        geom_point(data = subset(df_umap, marker_group == "None"), aes(color = marker_group), size = 0.6, alpha = 0.3) +
        geom_point(data = subset(df_umap, marker_group != "None"), aes(color = marker_group), size = 1, alpha = 0.9) +
        scale_color_manual(values = c("None" = "lightgrey", "Calb1" = "deepskyblue3", "Fos" = "deeppink2")) +
        labs(title = "Fos and Calb1 expression on UMAP (cluster-thresholded)") +
        plot_theme())
dev.off()

# 5) Neurotransmitter and Neuropeptide dotplots
p_nt <- DotPlot(subneurons, features = list(
  Glutamatergic = c("Slc17a6","Slc17a7"),
  GABAergic    = c("Gad1","Gad2"),
  Dopaminergic = c("Th")
), group.by = "seurat_clusters", dot.scale = 8) + plot_theme()
write_plot(p_nt, "06_dotplot_neurotransmitters", pdf=TRUE)

neuropep_list <- readLines("config/neuropeptide_genes.txt")
neuropep_list <- unique(neuropep_list[neuropep_list != ""])
p_np <- DotPlot(subneurons, features = list(Neuropeptides = neuropep_list), group.by = "seurat_clusters", dot.scale = 8) +
  plot_theme()
write_plot(p_np, "06_dotplot_neuropeptides", pdf=TRUE)

# 6) Top markers per cluster (3 genes) dotplot
Idents(subneurons) <- "seurat_clusters"

# filter out LOC and low pct.1, then take top 3 per cluster by p_adj then avg_log2FC
top3_table <- markers %>%
  filter(!grepl("^LOC", gene)) %>%
  filter(pct.1 >= 0.1) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

readr::write_csv(markers, "data/processed/06_subneurons_markers_all.csv")
readr::write_csv(top3_table, "data/processed/06_subneurons_markers_top3.csv")

gene_levels <- rev(unique(top3_table$gene))
tiff("figures/06_dotplot_top3_markers.tiff", width = 3000, height = 3000, res = 300)
print(DotPlot(subneurons, features = gene_levels, dot.scale = 8) + RotatedAxis() + coord_flip() + plot_theme())
dev.off()

saveRDS(subneurons, "data/processed/06_subneurons.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_06.txt")
