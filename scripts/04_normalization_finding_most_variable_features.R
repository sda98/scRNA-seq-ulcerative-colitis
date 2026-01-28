#!/usr/bin/env Rscript

# =============================================================================
# 04_normalization_finding_most_variable_features.R
#
# Purpose:
#   6) Normalize RNA counts using Seurat LogNormalize (scale.factor = 10,000)
#   7) Identify highly variable genes (HVGs) using vst (nfeatures = 2000)
#      - Print top 10 HVGs
#      - Save a styled HVG scatter plot (mean vs standardized variance)
#   8) Save the normalized + HVG-annotated Seurat object for downstream steps
#
# Inputs:
#   results/01_quality_control/objects/seurat_singlets.rds
#
# Outputs:
#   results/01_quality_control/plots/variable_features.png
#   results/01_quality_control/objects/seurat_norm_hvg.rds
#
# Notes:
#   - This script assumes doublets were removed in Script 03 (singlets object).
# =============================================================================

# =========================
# Libraries (required)
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(grid)    # for unit()
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "01_quality_control", "objects", "seurat_singlets.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 03_doublet_detection_scDblFinder.R first.")
}

out_dir    <- file.path(root_dir, "results", "01_quality_control")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load object
# =========================
cat("\nLoading singlet-only Seurat object...\n")
seurat_obj_filtered <- readRDS(in_obj)

# Basic sanity
if (!"RNA" %in% names(seurat_obj_filtered@assays)) {
  stop("Assay 'RNA' not found in the Seurat object. Check prior scripts.")
}

# =============================================================================
# 6) NORMALIZATION (LogNormalize)
# =============================================================================
cat("\n=== NORMALIZATION ===\n")
cat("Using LogNormalize (standard Seurat workflow)\n\n")
DefaultAssay(seurat_obj_filtered) <- "RNA"
seurat_obj_filtered <- NormalizeData(
  seurat_obj_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
cat("✓ Normalization complete\n")

# =============================================================================
# 7) FIND VARIABLE FEATURES
# =============================================================================
cat("\nFinding variable features...\n")
seurat_obj_filtered <- FindVariableFeatures(
  seurat_obj_filtered,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
cat("Top 10 variable features:\n")
top10 <- head(VariableFeatures(seurat_obj_filtered), 10)
print(top10)

# =============================================================================
# 7b) HVG PLOT (STYLED)
# =============================================================================
p0 <- VariableFeaturePlot(seurat_obj_filtered)
df_hvg <- p0$data
df_hvg$gene   <- rownames(df_hvg)
df_hvg$colors <- as.character(df_hvg$colors)
n_genes_total <- nrow(df_hvg)
other_label <- paste0("Least variable genes (", n_genes_total - 2000, ")")

p_hvg <- ggplot(df_hvg, aes(x = mean, y = variance.standardized)) +
  ggrepel::geom_label_repel(
    data = df_hvg[df_hvg$gene %in% top10, , drop = FALSE],
    aes(label = gene),
    fontface = "bold.italic",
    size = 9,
    fill = "#D6ECFF",
    alpha = 1,
    label.size = 0.6,
    label.r = unit(0.15, "lines"),
    label.padding = unit(0.25, "lines"),
    box.padding   = 1.9,
    point.padding = 0,
    segment.color = "black",
    segment.size = 2,
    segment.alpha = 1,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  geom_point(
    aes(fill = colors),
    shape = 21, color = "black", size = 6,
    stroke = 0.35, alpha = 1
  ) +
  scale_fill_manual(
    values = c("no" = "grey80", "yes" = "red",
               "FALSE" = "grey80", "TRUE" = "red"),
    breaks = c("no", "yes"),
    labels = c(other_label, "Most variable genes (2000)"),
    name = NULL,
    guide = guide_legend(
      override.aes = list(shape = 21, color = "black", size = 5, alpha = 1),
      byrow = TRUE
    )
  ) +
  guides(fill = guide_legend(title = NULL)) +
  scale_x_log10() +
  labs(
    x = "Average expression",
    y = "Standardized variance"
  ) +
  theme_classic(base_size = 50) +
  theme(
    axis.title = element_text(face = "bold", size = 55),
    axis.text = element_text(size = 55),
    legend.title = element_blank(),
    legend.text  = element_text(size = 35, face = "bold"),
    legend.spacing.y = unit(80, "pt"),
    legend.key.height = unit(80, "pt")
  )

ggsave(
  file.path(plot_dir, "variable_features.png"),
  p_hvg, width = 24, height = 18, dpi = 600
)

cat("✓ Variable features plot saved\n")

# =============================================================================
# 8) SAVE OBJECT FOR DOWNSTREAM SCRIPTS
# =============================================================================
out_obj <- file.path(object_dir, "seurat_norm_hvg.rds")
saveRDS(seurat_obj_filtered, out_obj)

cat("\n✓ Saved normalized + HVG object: ", out_obj, "\n")
cat("\n✓✓✓ NORMALIZATION + HVG COMPLETE ✓✓✓\n\n")

