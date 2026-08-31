# =============================================================================
# Fos_Calb1_quantification.R
#
# Fos / Calb1 expression across the final neuronal subclusters (subneurons.rds,
# stage 10 of PIL_snRNAseq_pipeline.R) and across the three experimental
# samples.
#
# EXPR_THRESHOLD and MIN_CELLS_PER_CLUSTER below are carried over from an
# earlier version of this analysis (log-normalized expression > 0.25, in at
# least 8 cells, for a cluster to be considered to reliably express the gene).
# CONFIRM these two numbers against the current manuscript Methods before
# treating this script's output as final -- they have not been verified
# against this manuscript revision in this session.
#
# Unlike that earlier version, a cluster failing the >=8-cell threshold is
# only EXCLUDED FROM THE "detected" CALL for that gene -- the underlying
# expression values in the Seurat object are never modified. Silently
# zeroing real measurements is not reproducible from the object alone and
# is difficult to justify under review.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("utils/plotting.R")   # defines plot_theme() and write_plot(plot, name, pdf = TRUE)
source("config/paths.R")     # defines out_dir

EXPR_THRESHOLD        <- 0.25   # log-normalized expression cutoff for "detected"
MIN_CELLS_PER_CLUSTER <- 8      # a cluster's call for a gene requires at least
                                # this many detected cells to be considered reliable

subneurons <- readRDS(file.path(out_dir, "subneurons.rds"))
DefaultAssay(subneurons) <- "RNA"
Idents(subneurons) <- "seurat_clusters"

genes_fc <- intersect(c("Fos", "Calb1"), rownames(subneurons))
stopifnot("Fos and/or Calb1 not found in the RNA assay" = length(genes_fc) == 2)

expr_data <- FetchData(subneurons, vars = c(genes_fc, "seurat_clusters", "orig.ident"))
expr_data$cell <- rownames(expr_data)

# --- per-cluster, per-gene reliability call (no values are modified) ---------
cluster_call <- expr_data %>%
  pivot_longer(all_of(genes_fc), names_to = "gene", values_to = "expr") %>%
  group_by(seurat_clusters, gene) %>%
  summarise(n_detected = sum(expr > EXPR_THRESHOLD),
            n_total    = n(),
            reliable   = n_detected >= MIN_CELLS_PER_CLUSTER,
            .groups = "drop")

write.csv(cluster_call, file.path(out_dir, "Fos_Calb1_cluster_reliability.csv"), row.names = FALSE)
message("--- CHECK: per-cluster detection counts for Fos/Calb1 ---")
print(cluster_call, n = Inf)

# --- per-sample proportions, within clusters called reliable -----------------
reliable_clusters <- cluster_call %>% filter(reliable) %>% select(seurat_clusters, gene)

per_sample <- expr_data %>%
  pivot_longer(all_of(genes_fc), names_to = "gene", values_to = "expr") %>%
  inner_join(reliable_clusters, by = c("seurat_clusters", "gene")) %>%
  mutate(detected = expr > EXPR_THRESHOLD) %>%
  group_by(gene, orig.ident) %>%
  summarise(n_detected = sum(detected), n_total = n(),
            pct = 100 * n_detected / n_total, .groups = "drop")

write.csv(per_sample, file.path(out_dir, "Fos_Calb1_by_sample.csv"), row.names = FALSE)
message("--- CHECK: Fos/Calb1 detection by sample (reliable clusters only) ---")
print(per_sample)

# --- violin plot: real expression values, all clusters ------------------------
if (!dir.exists(file.path(out_dir, "figures"))) dir.create(file.path(out_dir, "figures"), recursive = TRUE)

p_vln <- VlnPlot(subneurons, features = genes_fc, pt.size = 0, flip = TRUE,
                 group.by = "seurat_clusters") & plot_theme()
write_plot(p_vln, "violin_Fos_Calb1", pdf = TRUE)

# --- UMAP highlight, using the reliability call (not zeroed expression) ------
status <- expr_data %>%
  pivot_longer(all_of(genes_fc), names_to = "gene", values_to = "expr") %>%
  inner_join(cluster_call %>% select(seurat_clusters, gene, reliable),
            by = c("seurat_clusters", "gene")) %>%
  mutate(detected = expr > EXPR_THRESHOLD & reliable) %>%
  select(cell, gene, detected) %>%
  pivot_wider(names_from = gene, values_from = detected) %>%
  mutate(marker_group = case_when(
    Fos & Calb1 ~ "Both",
    Fos         ~ "Fos",
    Calb1       ~ "Calb1",
    TRUE        ~ "None"
  ))

df_umap <- as.data.frame(Embeddings(subneurons, "umap"))
umap_cols <- colnames(df_umap)   # "UMAP_1"/"UMAP_2" or "umap_1"/"umap_2" depending on Seurat version
df_umap$cell <- rownames(df_umap)
df_umap <- dplyr::left_join(df_umap, status[, c("cell", "marker_group")], by = "cell")

pal <- c(None = "lightgrey", Calb1 = "deepskyblue3", Fos = "deeppink2", Both = "darkorchid3")
p_hl <- ggplot(df_umap, aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]])) +
  geom_point(data = subset(df_umap, marker_group == "None"),
             aes(color = marker_group), size = 0.6, alpha = 0.3) +
  geom_point(data = subset(df_umap, marker_group != "None"),
             aes(color = marker_group), size = 1, alpha = 0.9) +
  scale_color_manual(values = pal) +
  labs(title = sprintf("Fos / Calb1 (expr > %.2f, clusters with >= %d detected cells)",
                       EXPR_THRESHOLD, MIN_CELLS_PER_CLUSTER)) +
  plot_theme()

write_plot(p_hl, "umap_Fos_Calb1", pdf = TRUE)

message("\nDone. Two files need your review before this is final:\n",
        "  - confirm EXPR_THRESHOLD / MIN_CELLS_PER_CLUSTER against the manuscript\n",
        "  - Fos_Calb1_by_sample.csv: nuclei from the same pooled library are not\n",
        "    independent observations of that condition -- report these percentages\n",
        "    descriptively, do not run a statistical test between samples on them.")
