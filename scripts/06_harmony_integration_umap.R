#!/usr/bin/env Rscript

# =============================================================================
# 06_harmony_integration_umap.R
#
# Purpose:
#   1) Visualize batch effects (sample-level variation) BEFORE integration
#      - Run UMAP on PCA embeddings (dims = 1:ncomp)
#      - Plot UMAP colored by sample_label
#   2) Perform Harmony integration (batch correction) using sample_id
#   3) Run standard workflow on Harmony embeddings
#      - UMAP + neighbors + clustering across multiple resolutions
#      - Set a default clustering resolution (0.4)
#   4) Visualize batch effects AFTER integration
#      - Plot UMAP colored by sample_label
#   5) Save integrated Seurat object for downstream scripts
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_pca.rds
#
# Outputs:
#   results/02_clustering_analysis/plots/umap_BEFORE_integration_by_sample.png
#   results/02_clustering_analysis/plots/umap_AFTER_integration_by_sample.png
#   results/02_clustering_analysis/objects/seurat_harmony_umap.rds
#
# Notes:
#   - This script expects metadata columns:
#       * gsm_id     (used to create sample_label)
#       * sample_id  (used as Harmony batch variable)
#   - ncomp is derived (in script 05_preprocess_pca_horns_parallel.R) via Horn's parallel analysis.
# =============================================================================

# =========================
# Libraries 
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "02_clustering_analysis", "objects", "seurat_pca.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 05_preprocess_pca_horns_parallel.R first.")
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
cat("\nLoading PCA-processed Seurat object...\n")
seurat_obj_filtered <- readRDS(in_obj)

DefaultAssay(seurat_obj_filtered) <- "RNA"

# Basic sanity checks
if (!"pca" %in% Reductions(seurat_obj_filtered)) {
  stop("PCA reduction not found in object. Make sure Script 05 ran RunPCA().")
}
if (!"gsm_id" %in% colnames(seurat_obj_filtered@meta.data)) {
  stop("Metadata column 'gsm_id' not found. Needed to create sample_label.")
}
if (!"sample_id" %in% colnames(seurat_obj_filtered@meta.data)) {
  stop("Metadata column 'sample_id' not found. Needed for Harmony integration.")
}

# =============================================================================
# 1) VISUALIZE BATCH EFFECTS BY SAMPLE (BEFORE INTEGRATION)
# =============================================================================
cat("\n=== BATCH EFFECT ASSESSMENT (BEFORE) ===\n")
cat("Visualizing sample-level variation before integration\n")

# ============================================
# UMAP BEFORE Integration
# ============================================
cat("\nRunning UMAP on non-integrated data...\n")
DefaultAssay(seurat_obj_filtered) <- "RNA"
seurat_obj_filtered <- RunUMAP(seurat_obj_filtered, dims = 1:ncomp, verbose = TRUE)

# ============================================
# Color by SAMPLE (12 samples)
# ============================================

# Create readable sample labels
seurat_obj_filtered$sample_label <- case_when(
  seurat_obj_filtered$gsm_id == "GSM7307094" ~ "UC Self-Control 1",
  seurat_obj_filtered$gsm_id == "GSM7307095" ~ "UC Self-Control 2",
  seurat_obj_filtered$gsm_id == "GSM7307096" ~ "UC Self-Control 3",
  seurat_obj_filtered$gsm_id == "GSM7307097" ~ "UC Self-Control 4",
  seurat_obj_filtered$gsm_id == "GSM7307098" ~ "Ulcerative Colitis 1",
  seurat_obj_filtered$gsm_id == "GSM7307099" ~ "Ulcerative Colitis 2",
  seurat_obj_filtered$gsm_id == "GSM7307100" ~ "Ulcerative Colitis 3",
  seurat_obj_filtered$gsm_id == "GSM7307101" ~ "Ulcerative Colitis 4",
  seurat_obj_filtered$gsm_id == "GSM7307102" ~ "Healthy Control 5",
  seurat_obj_filtered$gsm_id == "GSM7307103" ~ "Healthy Control 6",
  seurat_obj_filtered$gsm_id == "GSM7307104" ~ "Healthy Control 7",
  seurat_obj_filtered$gsm_id == "GSM7307105" ~ "Healthy Control 8",
  TRUE ~ "Unknown"
)

cat("Sample distribution:\n")
print(table(seurat_obj_filtered$sample_label))

# UMAP colored by sample
p_umap_sample <- DimPlot(
  seurat_obj_filtered,
  reduction = "umap",
  group.by = "sample_label",
  label = FALSE
) +
  ggtitle("UMAP before integration\n(by sample)") +
  labs(color = "Sample", x = "UMAP 1", y = "UMAP 2") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 30) +
  theme(
    plot.title   = element_text(face = "bold", size = 36),
    axis.title   = element_text(face = "bold", size = 34),
    axis.text    = element_text(size = 32),
    legend.title = element_text(face = "bold", size = 25),
    legend.text  = element_text(face = "bold", size = 18),
    legend.position = "right",
    legend.spacing.y = grid::unit(0.5, "cm"),
    legend.key.height = grid::unit(1, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

ggsave(
  file.path(plot_dir, "umap_BEFORE_integration_by_sample.png"),
  p_umap_sample, width = 14, height = 10, dpi = 600
)

cat("✓ Sample UMAP (before integration) saved\n")

# =============================================================================
# 2) HARMONY INTEGRATION
# =============================================================================
cat("\n=== HARMONY INTEGRATION ===\n")
cat("Using Harmony for gentle batch correction\n\n")

# Run Harmony on PCA
seurat_integrated <- RunHarmony(
  seurat_obj_filtered,
  group.by.vars = "sample_id",
  plot_convergence = FALSE,
  verbose = TRUE
)

cat("✓ Harmony integration complete\n")
cat("\nIntegrated reduction name: harmony\n")

# Ensure sample_label carries over
if (!"sample_label" %in% colnames(seurat_integrated@meta.data)) {
  seurat_integrated$sample_label <- seurat_obj_filtered$sample_label
}

# =============================================================================
# 3) STANDARD WORKFLOW ON INTEGRATED DATA
# =============================================================================
cat("\n=== STANDARD WORKFLOW (INTEGRATED) ===\n")
cat("Running UMAP + neighbors + clustering on Harmony embeddings...\n")

# Run UMAP on Harmony embeddings
seurat_integrated <- RunUMAP(
  seurat_integrated,
  reduction = "harmony",
  dims = 1:ncomp,
  verbose = TRUE
)
cat("✓ UMAP computed\n")

# Find neighbors using Harmony embeddings
seurat_integrated <- FindNeighbors(
  seurat_integrated,
  reduction = "harmony",
  dims = 1:ncomp,
  verbose = TRUE
)
cat("✓ Neighbors computed\n")

# Test multiple resolutions
set.seed(5)
seurat_integrated <- FindClusters(
  seurat_integrated,
  resolution = c(0.05, 0.1, 0.3, 0.4, 0.5, 0.8, 1),
  verbose = TRUE
)

# Check how many clusters each gives
cat("\nClusters at different resolutions:\n")
cat("Resolution 0.05:", length(unique(seurat_integrated$RNA_snn_res.0.05)), "\n")
cat("Resolution 0.1:",  length(unique(seurat_integrated$RNA_snn_res.0.1)),  "\n")
cat("Resolution 0.3:",  length(unique(seurat_integrated$RNA_snn_res.0.3)),  "\n")
cat("Resolution 0.4:",  length(unique(seurat_integrated$RNA_snn_res.0.4)),  "\n")
cat("Resolution 0.5:",  length(unique(seurat_integrated$RNA_snn_res.0.5)),  "\n")
cat("Resolution 0.8:",  length(unique(seurat_integrated$RNA_snn_res.0.8)),  "\n")
cat("Resolution 1:",    length(unique(seurat_integrated$RNA_snn_res.1)),    "\n")

Idents(seurat_integrated) <- "RNA_snn_res.0.4"

cat("\n✓ Clustering complete\n")
cat("Using resolution 0.4 as default\n")
cat("Number of clusters:", length(unique(Idents(seurat_integrated))), "\n\n")

# =============================================================================
# 4) UMAP AFTER Integration - BY SAMPLE
# =============================================================================
cat("\n=== BATCH EFFECT ASSESSMENT (AFTER) ===\n")
cat("Sample distribution (integrated object):\n")
print(table(seurat_integrated$sample_label))

p_umap_after_sample <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "sample_label",
  label = FALSE
) +
  ggtitle("UMAP after integration\n(by sample)") +
  labs(color = "Sample", x = "UMAP 1", y = "UMAP 2") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 30) +
  theme(
    plot.title   = element_text(face = "bold", size = 36),
    axis.title   = element_text(face = "bold", size = 34),
    axis.text    = element_text(size = 32),
    legend.title = element_text(face = "bold", size = 25),
    legend.text  = element_text(face = "bold", size = 18),
    legend.position = "right",
    legend.spacing.y = grid::unit(0.5, "cm"),
    legend.key.height = grid::unit(1, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

ggsave(
  file.path(plot_dir, "umap_AFTER_integration_by_sample.png"),
  p_umap_after_sample, width = 14, height = 10, dpi = 600
)

cat("✓ UMAP by sample (after integration) saved\n")

# =============================================================================
# 5) SAVE OBJECT FOR DOWNSTREAM SCRIPTS
# =============================================================================
out_obj <- file.path(object_dir, "seurat_harmony_umap.rds")
saveRDS(seurat_integrated, out_obj)

cat("\n✓ Saved integrated object: ", out_obj, "\n")
cat("\n✓✓✓ HARMONY + UMAP + CLUSTERING COMPLETE ✓✓✓\n\n")
