# 06_subneurons_analysis.R
# - subset neurons
# - re-SCTransform + PCA + UMAP + clustering 
# - per-sample proportions
# - exclude clusters based on:
#      (a) sample representation: any condition <10% AND highest >80%
#      (b) marker genes of adjacent brain regions
# - re-cluster PIL candidates
# - Calb1/Fos analysis
# - neurotransmitter + neuropeptide dotplots
# - top markers per cluster (3 genes) dotplot

source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

# 1) Subset neurons
obj <- readRDS("data/processed/05_annotated_all.rds")
 
# Neuronal clusters identified in step 05: GLU_1-4, GABA, DOPA
neuronal_idents <- c("GLU_1", "GLU_2", "GLU_3", "GLU_4", "GABA", "DOPA")
Idents(obj) <- "celltype_detailed"
neurons <- subset(obj, idents = neuronal_idents)
 
# Re-normalize and cluster neurons
# 30 PCs selected to capture finer transcriptional distinctions
# between closely related neuronal populations (explaining ~70% of variance)
neurons <- SCTransform(neurons, vst.flavor = "v2",
                       vars.to.regress = "nCount_RNA", verbose = FALSE)
neurons <- RunPCA(neurons, npcs = 50, verbose = FALSE)
neurons <- RunUMAP(neurons, dims = 1:30, min.dist = 0.4, spread = 1,
                   verbose = FALSE)
neurons <- FindNeighbors(neurons, dims = 1:30, verbose = FALSE)
neurons <- FindClusters(neurons, resolution = 1.5)
# → 37 neuronal clusters (0-36)
 
message("Neuronal clusters: ", length(levels(Idents(neurons))))
 
p_neurons <- DimPlot(neurons, label = TRUE) +
  ggtitle("Neuronal subclustering (dims=1:30, res=1.5)") + plot_theme()
write_plot(p_neurons, "06_umap_neurons_all", pdf = TRUE)
 
saveRDS(neurons, "data/processed/06_neurons_all.rds")

# 2) Per-sample proportions per cluster
prop_table <- prop.table(
  table(neurons$orig.ident, Idents(neurons)),
  margin = 2
) * 100
 
prop_df <- as.data.frame(prop_table)
colnames(prop_df) <- c("sample", "cluster", "percent")
readr::write_csv(prop_df, "data/processed/06_cluster_sample_percent.csv")

# 3) Identify clusters to EXCLUDE based on sample representation
#    Rule: exclude if ANY condition <10% AND highest condition >80%

imbalanced_clusters <- prop_df %>%
  group_by(cluster) %>%
  summarise(
    min_pct = min(percent),
    max_pct = max(percent),
    .groups = "drop"
  ) %>%
  filter(min_pct < 10 & max_pct > 80) %>%
  pull(cluster) %>%
  as.character()
 
message("Imbalanced clusters excluded: ",
        paste(imbalanced_clusters, collapse = ", "))
 

# 4) Identify clusters to EXCLUDE based on adjacent region marker genes
# Marker genes for adjacent brain regions
region_markers <- list(
  SN_dopaminerg    = c("Th", "Slc6a3"),
  SN_GABAerg       = c("Pax5", "Pou6f2", "Zfpm2"),      # Mendelsohn et al. 2024
  Cholinerg        = c("Chat", "Slc5a7"),
  Hippocampal      = c("Satb2", "Neurod6", "Slc17a7", "Prox1"),
  MG               = c("Synpo2", "Cnksr3", "Prkcd", "Rorb", "Gbx2"),
  ZI_hypothalamus  = c("Cdh23", "Pax6", "Pde11a"),
  MGE_int          = c("Lhx6", "Arx", "Dlx1"),
  Oligodendrocyta  = c("Mog", "Mbp")
)
 
# Find markers for all clusters
DefaultAssay(neurons) <- "RNA"
neurons <- NormalizeData(neurons, verbose = FALSE)
 
all_markers <- FindAllMarkers(
  neurons,
  only.pos    = TRUE,
  min.pct     = 0.2,
  logfc.threshold = 0.25,
  test.use    = "wilcox",
  verbose     = FALSE
)
 
# Identify clusters expressing regional marker genes
marker_excluded <- list()
for (region in names(region_markers)) {
  genes <- region_markers[[region]]
  hits <- all_markers %>%
    filter(gene %in% genes) %>%
    pull(cluster) %>%
    as.character() %>%
    unique()
  if (length(hits) > 0) {
    marker_excluded[[region]] <- hits
    message(region, ": clusters ", paste(hits, collapse = ", "))
  }
}
 
all_marker_excluded <- unique(unlist(marker_excluded))


# 5) Define PIL candidate clusters
all_clusters <- as.character(levels(Idents(neurons)))
 
excluded_clusters <- unique(c(imbalanced_clusters, all_marker_excluded))
PIL_candidates <- setdiff(all_clusters, excluded_clusters)
 
message("Excluded clusters: ", paste(sort(excluded_clusters), collapse = ", "))
message("PIL candidates: ", paste(sort(PIL_candidates), collapse = ", "))
 
# Build region assignment table
region_assignment <- data.frame(
  Cluster = c(
    "32",
    "1", "13", "22", "33", "35",
    "7", "23", "26", "36",
    "17",
    "9",
    "6",
    "21", "34",
    "25",
    paste(sort(PIL_candidates), collapse = ", ")
  ),
  Region = c(
    "ZI_hypothalamus",
    "MG", "MG", "MG", "MG", "MG",
    "Hippocampal", "Hippocampal", "Hippocampal", "Hippocampal",
    "SN_dopaminerg",
    "SN_GABAerg",
    "Oligodendrocyta",
    "Cholinerg", "Cholinerg",
    "MGE_int",
    "PIL_candidate"
  )
)
 
readr::write_csv(region_assignment,
                 "data/processed/06_region_assignment_TableS3.csv")
 

# 6) Subset PIL candidates and re-cluster → subneurons (21 subclusters)
subneurons <- subset(neurons, idents = PIL_candidates)
 
message("PIL candidate neurons: ", ncol(subneurons))
 
# Re-normalize and re-cluster
subneurons <- SCTransform(subneurons, vst.flavor = "v2",
                          vars.to.regress = "nCount_RNA", verbose = FALSE)
subneurons <- RunPCA(subneurons, npcs = 50, verbose = FALSE)
subneurons <- RunUMAP(subneurons, dims = 1:30, verbose = FALSE)
subneurons <- FindNeighbors(subneurons, dims = 1:30, verbose = FALSE)
subneurons <- FindClusters(subneurons,
                           graph.name = "SCT_snn",
                           resolution = 0.7)
 
message("PIL neuronal subclusters: ", length(levels(Idents(subneurons))))
 
p_sub <- DimPlot(subneurons, label = TRUE) +
  ggtitle("PIL subneurons (dims=1:30, res=0.7)") + plot_theme()
write_plot(p_sub, "06_umap_subneurons_PIL", pdf = TRUE)
 
p_sub_sample <- DimPlot(subneurons, group.by = "orig.ident") +
  ggtitle("PIL subneurons by condition") + plot_theme()
write_plot(p_sub_sample, "06_umap_subneurons_by_sample", pdf = TRUE)
 

# 7) Calb1/Fos analysis
#    Threshold: minimum 10 nuclei per cluster expressing gene > 0.25
#    (clusters with fewer cells considered statistically insufficient)

DefaultAssay(subneurons) <- "SCT"
expr_threshold  <- 0.25
min_cells       <- 10  # unified threshold for both Calb1 and Fos
 
calb1_expr <- GetAssayData(subneurons, layer = "data")["Calb1", ]
fos_expr   <- GetAssayData(subneurons, layer = "data")["Fos", ]
clusters   <- as.character(Idents(subneurons))
 
# Valid clusters
calb1_valid <- names(which(
  tapply(calb1_expr > expr_threshold, clusters, sum) >= min_cells))
fos_valid   <- names(which(
  tapply(fos_expr   > expr_threshold, clusters, sum) >= min_cells))
 
message("Calb1 valid clusters: ", paste(sort(calb1_valid), collapse = ", "))
message("Fos valid clusters: ",   paste(sort(fos_valid),   collapse = ", "))
 
# Per cluster × sample detailed table
meta <- subneurons@meta.data
meta$cell      <- rownames(meta)
meta$calb1_pos <- calb1_expr[meta$cell] > expr_threshold
meta$fos_pos   <- fos_expr[meta$cell]   > expr_threshold
 
result_detailed <- meta %>%
  dplyr::group_by(seurat_clusters, orig.ident) %>%
  dplyr::summarise(
    cell_count   = n(),
    fos_pos_n    = sum(fos_pos),
    calb1_pos_n  = sum(calb1_pos),
    double_pos_n = sum(fos_pos & calb1_pos),
    .groups = "drop"
  ) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(
    cluster_total = sum(cell_count),
    cell_pct      = round(100 * cell_count / cluster_total, 2),
    # Apply valid cluster filter
    fos_pos_n    = ifelse(as.character(seurat_clusters) %in%
                            fos_valid, fos_pos_n, 0),
    calb1_pos_n  = ifelse(as.character(seurat_clusters) %in%
                            calb1_valid, calb1_pos_n, 0),
    double_pos_n = ifelse(as.character(seurat_clusters) %in% fos_valid &
                            as.character(seurat_clusters) %in% calb1_valid,
                          double_pos_n, 0),
    # Per-sample percentages
    fos_pct_of_sample    = round(100 * fos_pos_n / cell_count, 2),
    calb1_pct_of_sample  = round(100 * calb1_pos_n / cell_count, 2),
    double_pct_of_sample = round(100 * double_pos_n / cell_count, 2),
    # Cluster-level percentages
    fos_pct_of_cluster    = round(100 * fos_pos_n / cluster_total, 2),
    calb1_pct_of_cluster  = round(100 * calb1_pos_n / cluster_total, 2),
    double_pct_of_cluster = round(100 * double_pos_n / cluster_total, 2)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(as.numeric(as.character(seurat_clusters)), orig.ident) %>%
  dplyr::select(
    Neuronal_subcluster             = seurat_clusters,
    Sample_name                     = orig.ident,
    Cell_count                      = cell_count,
    Cluster_total                   = cluster_total,
    `Cell_percentage(%)`            = cell_pct,
    cFos_positive_cells             = fos_pos_n,
    `cFos_of_sample(%)`             = fos_pct_of_sample,
    `cFos_of_cluster(%)`            = fos_pct_of_cluster,
    Calb1_positive_cells            = calb1_pos_n,
    `Calb1_of_sample(%)`            = calb1_pct_of_sample,
    `Calb1_of_cluster(%)`           = calb1_pct_of_cluster,
    Double_positive_cells           = double_pos_n,
    `Double_positive_of_sample(%)`  = double_pct_of_sample,
    `Double_positive_of_cluster(%)` = double_pct_of_cluster
  )
 
write.xlsx(result_detailed,
           "data/processed/06_Calb1_Fos_detailed_TableS6.xlsx",
           rowNames = FALSE)
 
# UMAP highlight: Calb1 and c-Fos
df_umap <- FetchData(subneurons, vars = c("umap_1", "umap_2"))
df_umap$cell <- rownames(df_umap)
 
expr_df <- data.frame(
  cell   = names(calb1_expr),
  Calb1  = as.numeric(calb1_expr),
  Fos    = as.numeric(fos_expr),
  cluster = clusters
)
 
expr_df <- expr_df %>%
  mutate(
    calb1_valid_flag = cluster %in% calb1_valid,
    fos_valid_flag   = cluster %in% fos_valid,
    marker_group = case_when(
      fos_valid_flag   & Fos   > expr_threshold ~ "Fos",
      calb1_valid_flag & Calb1 > expr_threshold ~ "Calb1",
      TRUE ~ "None"
    )
  )
 
df_umap <- dplyr::left_join(df_umap,
                             expr_df[, c("cell", "marker_group")],
                             by = "cell")
 
tiff("figures/06_umap_Calb1_cFos.tiff",
     width = 2500, height = 2000, res = 350)
print(ggplot(df_umap, aes(x = umap_1, y = umap_2)) +
        geom_point(data = subset(df_umap, marker_group == "None"),
                   aes(color = marker_group), size = 0.6, alpha = 0.3) +
        geom_point(data = subset(df_umap, marker_group != "None"),
                   aes(color = marker_group), size = 1, alpha = 0.9) +
        scale_color_manual(values = c("None"  = "lightgrey",
                                      "Calb1" = "deepskyblue3",
                                      "Fos"   = "deeppink2")) +
        labs(title = "Calb1 and c-Fos expression",
             color = "marker_group") +
        theme_minimal() +
        theme(panel.grid = element_blank(),
              axis.line  = element_line(color = "black"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 14),
              axis.text  = element_text(size = 12),
              legend.title = element_text(size = 13),
              legend.text  = element_text(size = 12)))
dev.off()
 
# Violin plot
tiff("figures/06_violin_Calb1_cFos.tiff",
     width = 3000, height = 1800, res = 300)
print(VlnPlot(subneurons,
              features = c("Calb1", "Fos"),
              assay    = "SCT",
              pt.size  = 0.5,
              group.by = "seurat_clusters") &
        theme(axis.text.x = element_text(size = 10)))
dev.off()
 
# Barplot: Calb1+, Fos+, double positive per cluster
result_cluster <- result_detailed %>%
  dplyr::group_by(Neuronal_subcluster) %>%
  dplyr::summarise(
    n_total  = sum(Cell_count),
    pct_calb1     = round(100 * sum(Calb1_positive_cells) / n_total, 1),
    pct_fos       = round(100 * sum(cFos_positive_cells)  / n_total, 1),
    pct_calb1_fos = round(100 * sum(Double_positive_cells)/ n_total, 1),
    .groups = "drop"
  ) %>%
  filter(pct_calb1 > 0 | pct_fos > 0)
 
result_plot <- result_cluster %>%
  tidyr::pivot_longer(
    cols      = c(pct_calb1, pct_fos, pct_calb1_fos),
    names_to  = "marker",
    values_to = "percent"
  ) %>%
  dplyr::mutate(
    marker = dplyr::recode(marker,
                           pct_calb1     = "Calb1+",
                           pct_fos       = "Fos+",
                           pct_calb1_fos = "Calb1+/Fos+"),
    marker = factor(marker, levels = c("Calb1+", "Fos+", "Calb1+/Fos+")),
    Neuronal_subcluster = factor(Neuronal_subcluster,
                                 levels = as.character(
                                   sort(as.numeric(
                                     unique(Neuronal_subcluster)))))
  )
 
tiff("figures/06_barplot_Calb1_cFos.tiff",
     width = 2400, height = 1800, res = 300)
print(ggplot(result_plot,
             aes(x = factor(Neuronal_subcluster), y = percent, fill = marker)) +
        geom_bar(stat = "identity", position = "dodge",
                 width = 0.7, color = "black", linewidth = 0.3) +
        scale_fill_manual(values = c("Calb1+"       = "deepskyblue3",
                                     "Fos+"         = "deeppink2",
                                     "Calb1+/Fos+"  = "purple3")) +
        labs(x = "Cluster", y = "% of neurons in cluster",
             fill = "Marker",
             title = "Proportion of Calb1+ and Fos+ neurons per cluster") +
        theme_classic() +
        theme(axis.text.x  = element_text(size = 12),
              axis.text.y  = element_text(size = 12),
              axis.title   = element_text(size = 13),
              legend.text  = element_text(size = 11),
              legend.title = element_text(size = 12),
              plot.title   = element_text(size = 14, face = "bold",
                                          hjust = 0.5)))
dev.off()
 

# 8) Neurotransmitter dotplot

tiff("figures/06_dotplot_neurotransmitters.tiff",
     width = 2000, height = 2000, res = 300)
print(DotPlot(subneurons,
              features = list(
                Glutamatergic = c("Slc17a6", "Slc17a7"),
                GABAergic     = c("Gad1", "Gad2", "Slc32a1")
              ),
              group.by = "seurat_clusters",
              dot.scale = 8,
              dot.min   = 0.05) +
        scale_y_discrete(limits = rev(levels(Idents(subneurons)))) +
        RotatedAxis() +
        theme(axis.text.x  = element_text(size = 14),
              strip.text.x = element_text(size = 14)) +
        guides(size  = guide_legend(title = "Percent\nExpressed"),
               color = guide_colorbar(title = "Average\nExpression")))
dev.off()
 
# 9) Neuropeptide dotplot

neuropep_genes <- c("Adcyap1", "Calca", "Calcb", "Cck", "Crh", "Grp",
                    "Npy", "Nts", "Pdyn", "Penk", "Pnoc", "Tac1", "Vgf")
 
tiff("figures/06_dotplot_neuropeptides.tiff",
     width = 2200, height = 2000, res = 300)
print(DotPlot(subneurons,
              features = neuropep_genes,
              group.by = "seurat_clusters",
              dot.scale = 8) +
        scale_y_discrete(limits = rev(levels(Idents(subneurons)))) +
        RotatedAxis() +
        theme(axis.text.x  = element_text(size = 12),
              axis.text.y  = element_text(size = 12)) +
        guides(size  = guide_legend(title = "Percent\nExpressed"),
               color = guide_colorbar(title = "Average\nExpression")) +
        ggtitle("Neuropeptides"))
dev.off()
 

# 10) FindAllMarkers for subneurons 

Idents(subneurons) <- "seurat_clusters"
 
subneuron_markers <- FindAllMarkers(
  subneurons,
  only.pos        = TRUE,
  min.pct         = 0.2,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)
 
readr::write_csv(subneuron_markers,
                 "data/processed/06_subneurons_markers_all_TableS4.csv")
 
# Top 3 markers per cluster dotplot
top3 <- subneuron_markers %>%
  filter(!grepl("^LOC", gene), pct.1 >= 0.2) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()
 
tiff("figures/06_dotplot_top3_markers.tiff",
     width = 3000, height = 3000, res = 300)
print(DotPlot(subneurons,
              features = rev(unique(top3$gene)),
              dot.scale = 8) +
        RotatedAxis() + coord_flip())
dev.off()
 

# Save
saveRDS(subneurons, "data/processed/06_subneurons.rds")
writeLines(capture.output(sessionInfo()),
           "data/processed/sessionInfo_06.txt")
 
message("Script 06 complete. Subneurons saved to data/processed/06_subneurons.rds")
