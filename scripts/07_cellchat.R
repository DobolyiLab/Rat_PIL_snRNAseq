# =============================================================================
# 07_cellchat.R  -- secreted signaling among neuronal subclusters
#
# FIXES vs original:
#  * reads 06_subneurons.rds (the original read 06_neurons.rds, which is never
#    written -> the script could not run)
#  * CellChat receives LOG-NORMALIZED data, as it requires; the old object's
#    RNA data slot still held raw counts
#  * selectK / identifyCommunicationPatterns are run before netAnalysis_river
#    and netAnalysis_dot, which otherwise error out
#  * CellChatDB.mouse used for rat data is now an explicit, logged decision
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(CellChat); library(NMF); library(ggalluvial) })

obj <- readRDS("data/processed/06_subneurons.rds")
DefaultAssay(obj) <- "RNA"
obj <- join_layers_if_split(obj, "RNA")
obj <- ensure_rna_normalized(obj)
data.input <- get_expr(obj, "RNA", "data")

lev  <- levels(droplevels(factor(obj$seurat_clusters)))
meta <- data.frame(labels = factor(paste0("cluster_", obj$seurat_clusters),
                                   levels = paste0("cluster_", lev)),
                   row.names = colnames(obj))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# Rat data against the mouse database: symbols are largely 1:1, but a minority
# of ligand-receptor pairs will not map. State this in the Methods.
message("[07] using CellChatDB.mouse for rat data -- orthology is assumed, not verified.")
data(CellChatDB.mouse)
cellchat@DB <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0.1)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

write_plot(
  netAnalysis_signalingRole_scatter(
    cellchat,
    dot.size = c(2, 6),
    label.size = 5,
    dot.alpha = 0.5
  ) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    ),
  "07_cellchat_scatter",
  pdf = TRUE
)

# --- communication patterns: required before river / dot plots ---------------
for (pat in c("outgoing", "incoming")) {
  pdf(sprintf("figures/07_cellchat_selectK_%s.pdf", pat), width = 8, height = 4)
  print(selectK(cellchat, pattern = pat))
  dev.off()
}
k_out <- params$cellchat$k_outgoing
k_in  <- params$cellchat$k_incoming
message(sprintf("[07] using k = %d (outgoing) and k = %d (incoming); ",
                k_out, k_in),
        "confirm these against the selectK elbow plots and adjust params.yml.")

# The pattern heatmap is a side effect of identifyCommunicationPatterns(), so
# the call is wrapped in a device to capture it. It must run ONCE per direction:
# NMF is expensive and a repeat run can settle on a different factorisation.
tiff("figures/07_cellchat_outgoing_patterns.tiff", width = 10, height = 8,
     units = "in", res = 300, compression = "lzw")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = k_out)
dev.off()

pdf("figures/07_cellchat_outgoing_river-dot.pdf", width = 10, height = 8)
print(netAnalysis_river(cellchat, pattern = "outgoing"))
print(netAnalysis_dot(cellchat, pattern = "outgoing"))
dev.off()

tiff("figures/07_cellchat_incoming_patterns.tiff", width = 10, height = 8,
     units = "in", res = 300, compression = "lzw")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = k_in)
dev.off()

pdf("figures/07_cellchat_incoming_river-dot.pdf", width = 10, height = 8)
print(netAnalysis_river(cellchat, pattern = "incoming"))
print(netAnalysis_dot(cellchat, pattern = "incoming"))
dev.off()

# ============================================================
# Signaling role heatmaps
# ============================================================

# PDF
pdf("figures/07_cellchat_role-heatmaps.pdf", width = 10, height = 8)

print(
  netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "outgoing",
    height = 11
  )
)

print(
  netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "incoming",
    height = 11
  )
)

dev.off()

# TIFF - outgoing
tiff(
  "figures/07_cellchat_role-heatmap_outgoing.tiff",
  width = 10,
  height = 8,
  units = "in",
  res = 300,
  compression = "lzw"
)

print(
  netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "outgoing",
    height = 11
  )
)

dev.off()

# TIFF - incoming
tiff(
  "figures/07_cellchat_role-heatmap_incoming.tiff",
  width = 10,
  height = 8,
  units = "in",
  res = 300,
  compression = "lzw"
)

print(
  netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "incoming",
    height = 11
  )
)

dev.off()

# ============================================================
# OUTGOING communication patterns
# ============================================================

# TIFF - outgoing river
tiff(
  "figures/07_cellchat_outgoing_river.tiff",
  width = 10,
  height = 8,
  units = "in",
  res = 300,
  compression = "lzw"
)

print(
  netAnalysis_river(
    cellchat,
    pattern = "outgoing"
  )
)

dev.off()

# TIFF - outgoing dot
tiff(
  "figures/07_cellchat_outgoing_dot.tiff",
  width = 10,
  height = 8,
  units = "in",
  res = 300,
  compression = "lzw"
)

print(
  netAnalysis_dot(
    cellchat,
    pattern = "outgoing"
  )
)

dev.off()

# ============================================================
# INCOMING communication patterns
# ============================================================

# TIFF - incoming river
tiff(
  "figures/07_cellchat_incoming_river.tiff",
  width = 10,
  height = 8,
  units = "in",
  res = 300,
  compression = "lzw"
)

print(
  netAnalysis_river(
    cellchat,
    pattern = "incoming"
  )
)

dev.off()

# TIFF - incoming dot
tiff(
  "figures/07_cellchat_incoming_dot.tiff",
  width = 10,
  height = 8,
  units = "in",
  res = 300,
  compression = "lzw"
)

print(
  netAnalysis_dot(
    cellchat,
    pattern = "incoming"
  )
)

dev.off()

readr::write_csv(subsetCommunication(cellchat), "data/processed/07_cellchat_interactions.csv")
saveRDS(cellchat, "data/processed/07_cellchat_object.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_07.txt")
