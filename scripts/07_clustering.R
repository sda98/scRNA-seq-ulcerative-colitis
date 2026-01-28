#!/usr/bin/env Rscript

# =============================================================================
# 07_clustering.R
#
# Purpose:
#   1) Generate publication-style UMAP plots colored by cluster
#      - res = 0.4 (with clean repelled labels)
#      - res = 0.3 (no labels) + dashed circle highlight
#      - res = 0.4 (no labels) + dashed circle highlight
#   2) Save intermediate clustered Seurat object (before annotations)
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_harmony_umap.rds
#
# Outputs:
#   results/02_clustering_analysis/plots/umap_clusters_res_0.4_clean.png
#   results/02_clustering_analysis/plots/umap_clusters_res_0.3_no_labels.png
#   results/02_clustering_analysis/plots/umap_clusters_res_0.4_no_labels.png
#   results/02_clustering_analysis/objects/seurat_integrated_clustered.rds
#
# Notes:
#   - Assumes the object already contains:
#       * UMAP reduction ("umap")
#       * clustering metadata columns (RNA_snn_res.*)
#   - This script only formats and exports cluster visualizations + saves an
#     intermediate object for downstream annotation scripts.
# =============================================================================

# =========================
# Libraries (required)
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(grid)   # for unit()
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "02_clustering_analysis", "objects", "seurat_harmony_umap.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 06_harmony_integration_umap.R first.")
}

out_dir    <- file.path(root_dir, "results", "02_clustering_analysis")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# Load object
# =========================
cat("\nLoading Harmony-integrated Seurat object...\n")
seurat_integrated <- readRDS(in_obj)

if (!"umap" %in% Reductions(seurat_integrated)) {
  stop("UMAP reduction not found in object. Make sure Script 06 ran RunUMAP().")
}

cat("✓ Object loaded\n")

# ============================================
# 15. UMAP plot (by cluster)
# ============================================

cat("\n=== CLUSTER UMAP EXPORTS ===\n")
cat("Generating UMAP plots colored by clusters\n")

# Define cluster colors matching cell type colors
# Cluster 4 gets a distinct color (grey for low-quality)
cluster_colors <- c(
  "0"  = "#E41A1C",
  "1"  = "#377EB8",
  "2"  = "#4DAF4A",
  "3"  = "#984EA3",
  "4"  = "#808080",
  "5"  = "#FF7F00",
  "6"  = "#FFFF33",
  "7"  = "#A65628",
  "8"  = "#F781BF",
  "9"  = "#999999",
  "10" = "#66C2A5",
  "11" = "#FC8D62",
  "12" = "#8DA0CB",
  "13" = "#E78AC3",
  "14" = "#A6D854",
  "15" = "#FFD92F",
  "16" = "#E5C494",
  "17" = "#B3B3B3",
  "18" = "#1B9E77"
)

# Get UMAP embeddings
umap_emb_cluster <- Embeddings(seurat_integrated, reduction = "umap") %>%
  as.data.frame()
colnames(umap_emb_cluster)[1:2] <- c("umap_1", "umap_2")

umap_coords_cluster <- umap_emb_cluster %>%
  mutate(
    cluster = as.character(Idents(seurat_integrated))
  )

# Calculate label positions
label_positions_cluster <- umap_coords_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

# Add colors and create lighter fill
label_positions_cluster$color <- cluster_colors[label_positions_cluster$cluster]

label_positions_cluster$fill_color <- sapply(label_positions_cluster$color, function(col) {
  rgb_vals <- col2rgb(col)
  lighter <- (rgb_vals + 200) / 255
  lighter[lighter > 1] <- 1
  rgb(lighter[1], lighter[2], lighter[3], alpha = 0.55)
})

p_umap_clusters <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  cols = cluster_colors,
  pt.size = 0.1
) +
  geom_label_repel(
    data = label_positions_cluster,
    aes(x = umap_1, y = umap_2, label = cluster, fill = I(fill_color)),
    size = 5,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE,
    max.overlaps = 50,
    box.padding = 0.2,
    point.padding = 0,
    segment.color = "black",
    segment.size = 0.6,
    min.segment.length = 0,
    label.padding = unit(0.25, "lines"),
    label.r = unit(0.15, "lines"),
    label.size = 0.25
  ) +
  ggtitle("UMAP\n(colored by cluster, res=0.4)") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 30) +
  theme(
    plot.title = element_text(face = "bold", size = 36),
    axis.title = element_text(face = "bold", size = 34),
    axis.text = element_text(size = 28),
    legend.position = "none"
  )

ggsave(
  file.path(plot_dir, "umap_clusters_res_0.4_clean.png"),
  p_umap_clusters, width = 10, height = 8, dpi = 600
)

cat("✓ Cluster UMAP saved\n")

# ============================================
# Save intermediate object (before annotations)
# ============================================

saveRDS(
  seurat_integrated,
  file = file.path(object_dir, "seurat_integrated_clustered.rds")
)

cat("\n✓ Intermediate object saved (before annotations)\n")
cat("\n✓✓✓ INTEGRATION & CLUSTERING COMPLETE ✓✓✓\n\n")

# ============================================
# 15A. UMAP with res = 0.3 - colors matched to res = 0.4
# ============================================

cat("\nGenerating UMAP for res = 0.3...\n")

Idents(seurat_integrated) <- "RNA_snn_res.0.3"
cat("Using resolution 0.3\n")
cat("Number of clusters:", length(unique(Idents(seurat_integrated))), "\n\n")

cluster_colors_res03 <- c(
  "0"  = "#377EB8",
  "1"  = "#E41A1C",
  "2"  = "#4DAF4A",
  "3"  = "#CDAD00",
  "4"  = "#984EA3",
  "5"  = "#808080",
  "6"  = "#FF7F00",
  "7"  = "#F781BF",
  "8"  = "#999999",
  "9"  = "#66C2A5",
  "10" = "#FC8D62",
  "11" = "#8DA0CB",
  "12" = "#E78AC3",
  "13" = "#A6D854",
  "14" = "#FFD92F",
  "15" = "#E5C494",
  "16" = "#1B9E77",
  "17" = "#00CED1"
)

# Get UMAP embeddings for cluster 3 to draw circle
umap_emb_03 <- Embeddings(seurat_integrated, reduction = "umap") %>%
  as.data.frame()
colnames(umap_emb_03)[1:2] <- c("umap_1", "umap_2")

umap_coords_03 <- umap_emb_03 %>%
  mutate(cluster = as.character(Idents(seurat_integrated)))

# Calculate center and radius for cluster 3
cluster3_coords <- umap_coords_03 %>%
  filter(cluster == "3")

cluster3_center <- c(
  x = median(cluster3_coords$umap_1),
  y = median(cluster3_coords$umap_2)
)

cluster3_distances <- sqrt(
  (cluster3_coords$umap_1 - cluster3_center["x"])^2 +
    (cluster3_coords$umap_2 - cluster3_center["y"])^2
)
cluster3_radius <- quantile(cluster3_distances, 0.95) + 0.5

# Create circle data (will be reused for both plots)
colonocyte_circle <- data.frame(
  x = cluster3_center["x"] + cluster3_radius * cos(seq(0, 2*pi, length.out = 100)),
  y = cluster3_center["y"] + cluster3_radius * sin(seq(0, 2*pi, length.out = 100))
)

p_umap_res03 <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  cols = cluster_colors_res03,
  pt.size = 0.1,
  label = FALSE
) +
  geom_path(
    data = colonocyte_circle,
    aes(x = x, y = y),
    color = "black",
    linewidth = 1.5,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggtitle("UMAP\n(colored by cluster, res=0.3)") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 30) +
  theme(
    plot.title = element_text(face = "bold", size = 36),
    axis.title = element_text(face = "bold", size = 34),
    axis.text = element_text(size = 28),
    legend.position = "none"
  )

ggsave(
  file.path(plot_dir, "umap_clusters_res_0.3_no_labels.png"),
  p_umap_res03, width = 10, height = 8, dpi = 600
)

cat("✓ UMAP (res = 0.3) cluster plot saved\n")

# ============================================
# 15B. UMAP with res = 0.4 with no labels
# ============================================

cat("\nGenerating UMAP for res = 0.4...\n")

Idents(seurat_integrated) <- "RNA_snn_res.0.4"
cat("Using resolution 0.4\n")
cat("Number of clusters:", length(unique(Idents(seurat_integrated))), "\n\n")

cluster_colors_res04 <- c(
  "0"  = "#E41A1C",
  "1"  = "#377EB8",
  "2"  = "#4DAF4A",
  "3"  = "#984EA3",
  "4"  = "#808080",
  "5"  = "#FF7F00",
  "6"  = "#FFFF33",
  "7"  = "#A65628",
  "8"  = "#F781BF",
  "9"  = "#999999",
  "10" = "#66C2A5",
  "11" = "#FC8D62",
  "12" = "#8DA0CB",
  "13" = "#E78AC3",
  "14" = "#A6D854",
  "15" = "#FFD92F",
  "16" = "#E5C494",
  "17" = "#B3B3B3",
  "18" = "#1B9E77"
)

p_umap_res04 <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  cols = cluster_colors_res04,
  pt.size = 0.1,
  label = FALSE
) +
  geom_path(
    data = colonocyte_circle,
    aes(x = x, y = y),
    color = "black",
    linewidth = 1.5,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggtitle("UMAP\n(colored by cluster, res=0.4)") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 30) +
  theme(
    plot.title = element_text(face = "bold", size = 36),
    axis.title = element_text(face = "bold", size = 34),
    axis.text = element_text(size = 28),
    legend.position = "none"
  )

ggsave(
  file.path(plot_dir, "umap_clusters_res_0.4_no_labels.png"),
  p_umap_res04, width = 10, height = 8, dpi = 600
)

cat("✓ UMAP (res = 0.4) cluster plot saved\n")
cat("\n✓✓✓ CLUSTER PLOTTING COMPLETE ✓✓✓\n\n")

