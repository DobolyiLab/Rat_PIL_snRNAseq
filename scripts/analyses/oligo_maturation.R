# =============================================================================
# analyses/oligo_maturation.R  -- is Oligo_NFOL actually a distinct maturation state?
#
# The markers in 05a are one-vs-ALL-CELLS: neurons and astrocytes sit in the
# denominator, so `Mbp` scoring high in an oligodendrocyte cluster only says
# "this is an oligodendrocyte", not where it sits in the maturation sequence.
# Restricting the comparison to the oligodendrocyte lineage isolates the
# maturation axis:  OPC -> COP -> NFOL -> MFOL -> MOL
# =============================================================================
source("scripts/00_setup.R")

obj <- readRDS("data/processed/05_annotated_all.rds")
DefaultAssay(obj) <- "RNA"
obj <- join_layers_if_split(obj, "RNA")
obj <- ensure_rna_normalized(obj)

ident_col <- if ("celltype_detailed" %in% colnames(obj@meta.data)) "celltype_detailed" else "seurat_clusters"
Idents(obj) <- ident_col

oligo_lineage <- grep("^Oligo|^OPC", levels(Idents(obj)), value = TRUE)
stopifnot("no oligodendrocyte-lineage clusters found" = length(oligo_lineage) > 1)
message("[oligo_maturation] lineage clusters: ", paste(oligo_lineage, collapse = ", "))

oli <- subset(obj, idents = oligo_lineage)
Idents(oli) <- droplevels(factor(Idents(oli)))
log_step("05c", oli, sprintf("(%d lineage clusters)", nlevels(Idents(oli))))

# --- DE within the lineage only ----------------------------------------------
mk <- FindAllMarkers(oli, assay = "RNA", only.pos = TRUE,
                     min.pct = 0.1, logfc.threshold = 0.25)
readr::write_csv(mk, "data/processed/qc/oligo_lineage_markers.csv")

top10 <- mk %>%
  filter(!grepl("^LOC|^AABR|^ENSRNOG", gene)) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
print(as.data.frame(top10))

# --- explicit maturation panel, in developmental order -----------------------
# The discriminating evidence for NFOL is the COMBINATION of early markers
# present (Enpp6, Anln) and mature markers absent (Opalin, Mal, Klk6).
stage_panel <- list(
  OPC   = c("Pdgfra", "Cspg4", "Vcan", "Ptprz1"),
  COP   = c("Gpr17", "Bmp4", "Neu4", "Sema5b"),
  NFOL  = c("Enpp6", "Anln", "Arhgef37", "Rasgrp3", "Tmem2", "Itpr2"),
  MFOL  = c("Mal", "Ctps1", "Opalin", "Mog"),
  MOL   = c("Klk6", "Apod", "Ptgds", "Il33", "Car2", "Hapln2", "Trf", "Tf")
)
stage_panel <- lapply(stage_panel, function(g) present_features(oli, g, "RNA", "stage panel"))
stage_panel <- stage_panel[lengths(stage_panel) > 0]

write_plot(DotPlot(oli, features = stage_panel, group.by = ident_col,
                   dot.scale = 8, dot.min = params$plots$dot_min) +
             plot_theme() + rotated_x(),
           "fig_oligo_maturation", pdf = TRUE)

# --- per-cluster mean expression of each stage panel -------------------------
# A cluster's stage is the panel it scores highest on. Reported as the mean
# fraction of the panel's genes detected per nucleus, so panels of different
# size stay comparable.
cts <- get_expr(oli, "RNA", "counts")
grp <- factor(oli@meta.data[[ident_col]])
stage_score <- vapply(stage_panel, function(g) {
  pos <- Matrix::colMeans(cts[g, , drop = FALSE] > 0)
  tapply(pos, grp, mean) * 100
}, numeric(nlevels(grp)))

stage_tbl <- as.data.frame(round(stage_score, 1))
stage_tbl$cluster <- rownames(stage_tbl)
stage_tbl$n_cells <- as.integer(table(grp))
stage_tbl$best_stage <- colnames(stage_score)[apply(stage_score, 1, which.max)]
stage_tbl <- stage_tbl %>% select(cluster, n_cells, best_stage, everything())

readr::write_csv(stage_tbl, "data/processed/qc/oligo_stage_scores.csv")
print(as.data.frame(stage_tbl))

message("\n[oligo_maturation] A cluster only earns a stage name when the panel it scores highest on ",
        "also matches the ABSENCE of the later markers. If two clusters share a ",
        "best_stage, name them neutrally (Oligo_1, Oligo_2) and describe the ",
        "markers in the text instead.")

writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_oligo_maturation.txt")
