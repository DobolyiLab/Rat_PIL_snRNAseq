# Step 05 - Cell type annotation for all cells
# This script generates UMAP and dot plots for all cells and optionally renames clusters
# based on a CSV mapping file.

source("scripts/00_setup.R")

# Load the processed single-cell object after SoupX and SCTransform
obj <- readRDS("data/processed/04_singlets_soupx_sct.rds")

# Set RNA assay as default and normalize if data slot is empty
DefaultAssay(obj) <- "RNA"
if (max(GetAssayData(obj, assay = "RNA", slot = "data")) <= 0) {
  obj <- NormalizeData(obj, assay = "RNA")
}

# Load marker genes used in figures
markers <- readr::read_csv("config/markers_used_in_figures.csv", show_col_types = FALSE)$gene

# Plot dotplot for all markers across all clusters
write_plot(
  DotPlot(obj, features = markers, cols = c("lightgrey", "midnightblue"),
          col.min = -1, scale.min = 0) + plot_theme(),
  "05_dotplot_markers_all", pdf = TRUE
)

# Plot UMAPs by cluster and by sample
write_plot(
  DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + plot_theme(),
  "05_umap_by_cluster_all", pdf = TRUE
)
write_plot(
  DimPlot(obj, reduction = "umap", group.by = "orig.ident", label = TRUE) + plot_theme(),
  "05_umap_by_sample_all", pdf = TRUE
)

# Cluster renaming based on CSV mapping file
cur <- sort(as.character(unique(obj$seurat_clusters)))
if (!file.exists("config/cluster_rename.csv")) {
  # If mapping file does not exist, create a template
  readr::write_csv(data.frame(cluster_id = cur, new_name = ""), "config/cluster_rename_template.csv")
  message("Fill in 'config/cluster_rename_template.csv' and save it as 'cluster_rename.csv', then re-run this script.")
} else {
  # Read the mapping file and check columns
  ren <- readr::read_csv("config/cluster_rename.csv", show_col_types = FALSE)
  stopifnot(all(c("cluster_id", "new_name") %in% names(ren)))
  
  # Keep only valid mappings
  ren <- ren[ren$cluster_id %in% cur & ren$new_name != "", , drop = FALSE]
  
  # Apply new names
  Idents(obj) <- "seurat_clusters"
  map <- setNames(ren$new_name, ren$cluster_id)
  obj$celltype_detailed <- factor(ifelse(
    as.character(Idents(obj)) %in% names(map),
    map[as.character(Idents(obj))],
    as.character(Idents(obj))
  ))
  Idents(obj) <- "celltype_detailed"
  
  # Plots after renaming
  write_plot(
    DimPlot(obj, reduction = "umap", label = TRUE) + plot_theme(),
    "05_umap_celltype_labeled_all", pdf = TRUE
  )
  write_plot(
    DotPlot(obj, features = markers, group.by = "celltype_detailed",
            cols = c("lightgrey", "midnightblue"), col.min = -1, scale.min = 0) + plot_theme(),
    "05_dotplot_by_celltype_all", pdf = TRUE
  )
  
  # Save annotated object
  saveRDS(obj, "data/processed/05_annotated_all.rds")
}
