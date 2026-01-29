#!/usr/bin/env Rscript

# =============================================================================
# 13_final_cell_annotation.R
#
# Purpose:
#   1) Load annotated Seurat object with manual cell type labels.
#   2) Remove the low-quality population labeled (Cluster 4).
#   3) Add broad cell type categories (celltype_broad).
#   4) Generate publication-style UMAPs:
#        - Final cell type labels
#        - Broad categories
#        - Final cell types split by condition (with colonocyte highlight)
#   5) Save the cleaned Seurat object for downstream scripts.
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_annotated_singler.rds
#     (or change `in_obj` below to whatever annotated object you want to use)
#
# Outputs:
#   results/02_clustering_analysis/sessionInfo_13.txt
#   results/02_clustering_analysis/plots/umap_celltype_final.png
#   results/02_clustering_analysis/plots/umap_celltype_broad.png
#   results/02_clustering_analysis/plots/umap_celltype_by_condition.png
#   results/02_clustering_analysis/objects/seurat_integrated_clean.rds
#
# Notes:
#   - Broad categories are mapped from `celltype_final` (so names must match).
#   - Condition-split UMAP requires a metadata column named `condition`.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggh4x)
  library(ggforce)
  library(grid)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(
  root_dir, "results", "02_clustering_analysis", "objects",
  "seurat_integrated_annotated_singler.rds"
)
if (!file.exists(in_obj)) {
  stop(
    "Input object not found: ", in_obj,
    "\nUpdate `in_obj` to the correct annotated object (e.g., seurat_integrated_annotated.rds)."
  )
}

out_dir    <- file.path(root_dir, "results", "02_clustering_analysis")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_13.txt"))

# =========================
# 1) LOAD SEURAT OBJECT
# =========================
cat("\nLoading annotated Seurat object...\n")
seurat_integrated <- readRDS(in_obj)

if (!"celltype_final" %in% colnames(seurat_integrated@meta.data)) {
  stop("Missing 'celltype_final' in metadata. This script expects manual annotations already added.")
}
if (!"umap" %in% Reductions(seurat_integrated)) {
  stop("UMAP reduction not found. Make sure UMAP exists in this object.")
}

cat("✓ Object loaded\n")

# ============================================
# 2) REMOVE CLUSTER 4 (Low quality, see script 11)
# ============================================
cat("\n=== Removing low-quality cluster ===\n")

cells_before <- ncol(seurat_integrated)

seurat_clean <- subset(
  seurat_integrated,
  subset = celltype_final != "Low-quality (excluded)"
)

cells_after <- ncol(seurat_clean)

cat("Cells before:", cells_before, "\n")
cat("Cells after:",  cells_after,  "\n")
cat("Removed:",      cells_before - cells_after, "low-quality cells\n")

# ============================================
# 3) ADDING BROAD CATEGORIES
# ============================================
cat("\nAdding broad cell type categories...\n")

# Rename colonocyte labels to preferred naming
seurat_clean$celltype_final <- gsub(
  "Colonocytes (less differentiated)", 
  "Colonocytes (less mature)", 
  seurat_clean$celltype_final,
  fixed = TRUE
)
seurat_clean$celltype_final <- gsub(
  "Colonocytes (more differentiated)", 
  "Colonocytes (more mature)", 
  seurat_clean$celltype_final,
  fixed = TRUE
)

broad_categories <- c(
  "CD4+ T cells"                               = "Immune - T cells",
  "Cytotoxic lymphocytes (CD8+ T / NK)"        = "Immune - T cells",
  "B cells (naive/memory)"                     = "Immune - B/Plasma",
  "Germinal center B cells"                    = "Immune - B/Plasma",
  "Cycling B cells"                            = "Immune - B/Plasma",
  "Plasma cells (IgA+)"                        = "Immune - B/Plasma",
  "Plasma cells (IgG+)"                        = "Immune - B/Plasma",
  "Macrophages"                                = "Immune - Myeloid",
  "Mast cells"                                 = "Immune - Myeloid",
  "Colonocytes (less mature)"                  = "Epithelial",
  "Colonocytes (more mature)"                  = "Epithelial",
  "Goblet cells"                               = "Epithelial",
  "Tuft cells"                                 = "Epithelial",
  "Enteroendocrine cells (L-cell enriched)"    = "Epithelial",
  "Inflammatory fibroblasts"                   = "Stromal",
  "Endothelial cells"                          = "Stromal",
  "Pericytes / SMCs"                           = "Stromal",
  "Enteric glia"                               = "Neural"
)

# Create named vector with cell barcodes
broad_labels <- broad_categories[as.character(seurat_clean$celltype_final)]
names(broad_labels) <- colnames(seurat_clean)

seurat_clean <- AddMetaData(
  seurat_clean,
  metadata = broad_labels,
  col.name = "celltype_broad"
)

seurat_clean$celltype_broad <- factor(
  seurat_clean$celltype_broad,
  levels = c("Immune - T cells", "Immune - B/Plasma", "Immune - Myeloid",
             "Epithelial", "Stromal", "Neural")
)

cat("✓ Broad categories added\n")
print(table(seurat_clean$celltype_broad, useNA = "ifany"))

# Adding colors

celltype_colors <- c(
  "CD4+ T cells"                               = "#E41A1C",
  "Plasma cells (IgA+)"                        = "#377EB8",
  "Cytotoxic lymphocytes (CD8+ T / NK)"        = "#4DAF4A",
  "B cells (naive/memory)"                     = "#984EA3",
  "Plasma cells (IgG+)"                        = "#FF7F00",
  "Colonocytes (less mature)"                  = "#FFFF33",
  "Colonocytes (more mature)"                  = "#A65628",
  "Inflammatory fibroblasts"                   = "#F781BF",
  "Macrophages"                                = "#999999",
  "Cycling B cells"                            = "#66C2A5",
  "Goblet cells"                               = "#FC8D62",
  "Endothelial cells"                          = "#8DA0CB",
  "Mast cells"                                 = "#E78AC3",
  "Enteric glia"                               = "#A6D854",
  "Tuft cells"                                 = "#FFD92F",
  "Pericytes / SMCs"                           = "#E5C494",
  "Germinal center B cells"                    = "#B3B3B3",
  "Enteroendocrine cells (L-cell enriched)"    = "#1B9E77"
)

broad_colors <- c(
  "Immune - T cells"    = "#E41A1C",
  "Immune - B/Plasma"   = "#377EB8",
  "Immune - Myeloid"    = "#4DAF4A",
  "Epithelial"          = "#984EA3",
  "Stromal"             = "#FF7F00",
  "Neural"              = "#A6D854"
)

# Helper: lighter fill colors for labels
.light_fill <- function(col_hex) {
  rgb_vals <- grDevices::col2rgb(col_hex)
  lighter <- (rgb_vals + 200) / 255
  lighter[lighter > 1] <- 1
  grDevices::rgb(lighter[1], lighter[2], lighter[3], alpha = 0.7)
}

# ============================================
# 4) UMAP - FINAL CELL TYPES
# ============================================
cat("\nGenerating final cell type UMAP...\n")

umap_emb <- Embeddings(seurat_clean, reduction = "umap") %>% as.data.frame()
colnames(umap_emb)[1:2] <- c("umap_1", "umap_2")

umap_coords <- umap_emb %>%
  mutate(celltype = seurat_clean$celltype_final)

label_positions <- umap_coords %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarize(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

label_positions$color <- celltype_colors[as.character(label_positions$celltype)]
label_positions$fill_color <- vapply(label_positions$color, .light_fill, character(1))
label_positions$celltype_wrapped <- gsub("\\s*\\(", "\n(", label_positions$celltype)

p_umap_final <- DimPlot(
  seurat_clean,
  reduction = "umap",
  group.by  = "celltype_final",
  cols      = celltype_colors,
  pt.size   = 0.1
) +
  geom_label_repel(
    data = label_positions,
    aes(x = umap_1, y = umap_2, label = celltype_wrapped, fill = I(fill_color)),
    size = 3,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE,
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0,
    nudge_y = 0,
    segment.color = "black",
    segment.size = 0.4,
    min.segment.length = 0,
    label.padding = unit(0.25, "lines"),
    label.r = unit(0.15, "lines"),
    label.size = 0.25
  ) +
  ggtitle("Cell Type Annotation") +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
  ) +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 30) +
  theme(
    plot.title   = element_text(face = "bold", size = 33),
    axis.title   = element_text(face = "bold", size = 31),
    axis.text    = element_text(size = 28),
    plot.caption = element_text(size = 14, hjust = 0, face = "italic"),
    legend.position = "none"
  )

ggsave(
  file.path(plot_dir, "umap_celltype_final.png"),
  p_umap_final, width = 12, height = 10, dpi = 600, bg = "white"
)
cat("✓ Final annotation UMAP saved\n")

# ============================================
# UMAP - Broad categories
# ============================================
cat("\nGenerating broad category UMAP...\n")

umap_coords_broad <- umap_emb %>%
  mutate(celltype_broad = seurat_clean$celltype_broad)

label_positions_broad <- umap_coords_broad %>%
  dplyr::group_by(celltype_broad) %>%
  dplyr::summarize(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

label_positions_broad$color <- broad_colors[as.character(label_positions_broad$celltype_broad)]
label_positions_broad$fill_color <- vapply(label_positions_broad$color, .light_fill, character(1))

p_umap_broad <- DimPlot(
  seurat_clean,
  reduction = "umap",
  group.by  = "celltype_broad",
  cols      = broad_colors,
  pt.size   = 0.1
) +
  geom_label_repel(
    data = label_positions_broad,
    aes(x = umap_1, y = umap_2, label = celltype_broad, fill = I(fill_color)),
    size = 5,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0,
    segment.color = "black",
    segment.size = 0.7,
    min.segment.length = 0,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.15, "lines"),
    label.size = 0.3
  ) +
  ggtitle("Broad Cell Type Categories") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 30) +
  theme(
    plot.title   = element_text(face = "bold", size = 33),
    axis.title   = element_text(face = "bold", size = 31),
    axis.text    = element_text(size = 28),
    legend.position = "none"
  )

ggsave(
  file.path(plot_dir, "umap_celltype_broad.png"),
  p_umap_broad, width = 12, height = 10, dpi = 600, bg = "white"
)
cat("✓ Broad category UMAP saved\n")

# ============================================
# UMAP - Split by condition (if present)
# ============================================
cat("\nGenerating condition-split UMAP...\n")

if (!"condition" %in% colnames(seurat_clean@meta.data)) {
  cat("⚠ 'condition' column not found in metadata. Skipping condition-split UMAP.\n")
} else {
  
  umap_coords_condition <- umap_emb %>%
    mutate(
      celltype   = seurat_clean$celltype_final,
      condition  = seurat_clean$condition
    )
  
  umap_coords_condition$condition <- factor(
    umap_coords_condition$condition,
    levels = c("Healthy Control", "UC Self-Control", "Ulcerative Colitis")
  )
  
  condition_labels <- c(
    "Healthy Control"     = "Healthy Control",
    "UC Self-Control"     = "UC (Self-Control)",
    "Ulcerative Colitis"  = "UC (Inflamed)"
  )
  
  # Dashed circle around colonocytes (clusters 6 and 7 label names)
  colonocyte_cells <- umap_coords_condition %>%
    filter(celltype %in% c("Colonocytes (less mature)", "Colonocytes (more mature)"))
  
  circle_center_x <- mean(colonocyte_cells$umap_1)
  circle_center_y <- mean(colonocyte_cells$umap_2)
  
  distances <- sqrt((colonocyte_cells$umap_1 - circle_center_x)^2 + (colonocyte_cells$umap_2 - circle_center_y)^2)
  circle_radius <- as.numeric(stats::quantile(distances, 0.95)) * 1.1
  
  cat("Circle center:", circle_center_x, ",", circle_center_y, "\n")
  cat("Circle radius:", circle_radius, "\n")
  cat("Number of colonocytes found:", nrow(colonocyte_cells), "\n")
  
  circle_data <- data.frame(
    x0 = circle_center_x,
    y0 = circle_center_y,
    r  = circle_radius,
    condition = factor(
      c("Healthy Control", "UC Self-Control", "Ulcerative Colitis"),
      levels = c("Healthy Control", "UC Self-Control", "Ulcerative Colitis")
    )
  )
  
  p_umap_condition <- ggplot(umap_coords_condition, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = celltype_colors) +
    ggforce::geom_circle(
      data = circle_data,
      aes(x0 = x0, y0 = y0, r = r),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 1.2,
      color = "black"
    ) +
    ggh4x::facet_wrap2(
      ~ condition,
      ncol = 3,
      labeller = labeller(condition = condition_labels),
      scales = "fixed",
      axes = "all"
    ) +
    labs(title = NULL, x = "UMAP 1", y = "UMAP 2") +
    scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
    scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 28) +
    theme(
      axis.title = element_text(face = "bold", size = 40),
      axis.title.x = element_text(margin = margin(t = 20, unit = "pt")),
      axis.title.y = element_text(margin = margin(r = 20, unit = "pt")),
      axis.text = element_text(size = 35),
      axis.ticks = element_line(linewidth = 1.2),
      axis.ticks.length = unit(0.3, "cm"),
      axis.line = element_line(linewidth = 1),
      strip.text = element_text(face = "bold", size = 32),
      strip.background = element_rect(fill = "grey90", color = "black"),
      legend.position = "none",
      panel.spacing = unit(1.5, "lines")
    )
  
  ggsave(
    file.path(plot_dir, "umap_celltype_by_condition.png"),
    p_umap_condition, width = 24, height = 10, dpi = 600, bg = "white"
  )
  cat("✓ Condition-split UMAP saved\n")
}

# ============================================
# 5) SAVE SEURAT OBJECT
# ============================================
out_obj <- file.path(object_dir, "seurat_integrated_clean.rds")
saveRDS(seurat_clean, out_obj)
cat("\n✓ Saved cleaned Seurat object: ", out_obj, "\n")
