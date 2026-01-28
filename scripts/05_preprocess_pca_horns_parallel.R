#!/usr/bin/env Rscript

# =============================================================================
# 05_preprocess_pca_horns_parallel.R
#
# Purpose:
#   1) Scale data (ScaleData) using previously selected HVGs
#   2) Run PCA (RunPCA) using HVGs
#   3) Estimate an "optimal" number of PCs using Horn’s Parallel Analysis
#      (sampling up to 2000 cells and 500 HVGs for computational efficiency)
#   4) Save a styled elbow plot with a vertical line at the Horn-derived PC count
#   5) Save the PCA-processed Seurat object for downstream clustering/integration
#
# Inputs:
#   results/01_quality_control/objects/seurat_norm_hvg.rds
#
# Outputs:
#   results/02_clustering_analysis/plots/qc_PCA_elbow.png
#   results/02_clustering_analysis/objects/seurat_pca.rds
#
# Notes:
#   - Horn’s parallel analysis is computed on a sampled subset of cells/genes to
#     keep runtime reasonable.
#   - The elbow plot is based on PCA standard deviations computed by Seurat.
# =============================================================================

# =========================
# Libraries (required)
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "01_quality_control", "objects", "seurat_norm_hvg.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 04_normalization_finding_most_variable_features.R first.")
}

# Input comes from QC stage
in_stage  <- file.path(root_dir, "results", "01_quality_control")

# Output goes to clustering stage
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
cat("\nLoading normalized + HVG Seurat object...\n")
seurat_obj_filtered <- readRDS(in_obj)

# Basic sanity
if (length(VariableFeatures(seurat_obj_filtered)) == 0) {
  stop("No VariableFeatures found. Make sure Script 04 ran FindVariableFeatures().")
}
DefaultAssay(seurat_obj_filtered) <- "RNA"

# =============================================================================
# 1) SCALING + PCA
# =============================================================================
cat("\n=== SCALING AND PCA ===\n")

seurat_obj_filtered <- ScaleData(
  seurat_obj_filtered,
  features = VariableFeatures(seurat_obj_filtered),
  verbose = FALSE
)
cat("✓ Data scaled\n")

seurat_obj_filtered <- RunPCA(
  seurat_obj_filtered,
  features = VariableFeatures(seurat_obj_filtered),
  verbose = FALSE
)
cat("✓ PCA computed\n")

# =============================================================================
# 2) Horn's Parallel Analysis to determine optimal PCs
# =============================================================================
cat("\nRunning Horn's Parallel Analysis...\n")

set.seed(42)
n_cells_sample <- min(2000, ncol(seurat_obj_filtered))
cells_use <- sample(colnames(seurat_obj_filtered), n_cells_sample)
genes_use <- head(VariableFeatures(seurat_obj_filtered), 500)

X <- GetAssayData(seurat_obj_filtered, assay = "RNA", layer = "scale.data")[genes_use, cells_use]
X <- t(as.matrix(X))  # cells x genes

# -----------------------------
# Remove problematic genes (NA/Inf or sd=0)
# -----------------------------
genes_before <- ncol(X)
bad_nf <- colSums(!is.finite(X)) > 0   # has NA/Inf
bad_sd <- apply(X, 2, sd) == 0         # zero variance
X <- X[, !(bad_nf | bad_sd), drop = FALSE]

cat("Genes before cleaning:", genes_before, "\n")
cat("Genes dropped:", sum(bad_nf | bad_sd), "\n")
cat("  - dropped (non-finite):", sum(bad_nf), "\n")
cat("  - dropped (sd=0):", sum(bad_sd), "\n")
cat("Genes remaining:", ncol(X), "\n\n")

# -----------------------------
# Horn parallel analysis
# -----------------------------
cat("Computing eigenvalues...\n")
obs <- eigen(cor(X), only.values = TRUE)$values

cat("Running", 30, "permutations...\n")
n_iter <- 30
sim <- replicate(n_iter, {
  R <- matrix(rnorm(nrow(X) * ncol(X)), nrow = nrow(X))
  eigen(cor(R), only.values = TRUE)$values
})

thr <- apply(sim, 1, quantile, probs = 0.95)
ncomp <- sum(obs > thr)

cat("✓ Optimal number of PCs:", ncomp, "\n\n")

# =============================================================================
# 3) Styled Elbow Plot with Horn ncomp line
# =============================================================================
cat("Creating styled elbow plot...\n")

p_elbow <- ElbowPlot(seurat_obj_filtered, ndims = 50)
y_top <- max(p_elbow$data$stdev, na.rm = TRUE) * 1.05
lab   <- as.expression(bquote(PC[optimal] == .(ncomp)))

p_elbow <- p_elbow +
  geom_vline(
    xintercept = ncomp,
    linetype   = "dashed",
    color      = "red",
    linewidth  = 2
  ) +
  annotate(
    "text",
    x = ncomp,
    y = y_top,
    label = lab,
    hjust = 1.1,
    vjust = 1,
    size  = 10,
    fontface = "bold",
    color = "black"
  ) +
  ggtitle("Elbow plot (PCA)") +
  labs(
    x = "Principal Components",
    y = "Standard Deviation"
  ) +
  theme_classic(base_size = 32) +
  theme(
    plot.title = element_text(face = "bold", size = 36),
    axis.title = element_text(face = "bold", size = 34),
    axis.text  = element_text(size = 30)
  )

ggsave(
  file.path(plot_dir, "qc_PCA_elbow.png"),
  p_elbow, width = 12, height = 10, dpi = 600
)

cat("✓ Elbow plot saved\n")
cat("\nRecommended PCs for downstream analysis:", ncomp, "\n")

# =============================================================================
# 4) SAVE OBJECT FOR DOWNSTREAM SCRIPTS
# =============================================================================
out_obj <- file.path(object_dir, "seurat_pca.rds")
saveRDS(seurat_obj_filtered, out_obj)

cat("\n✓ Saved PCA object: ", out_obj, "\n")
cat("\n✓✓✓ PCA + HORN PARALLEL COMPLETE ✓✓✓\n\n")

cat("\n✓ Saved PCA object: ", out_obj, "\n")
cat("\n✓✓✓ PCA + HORN PARALLEL COMPLETE ✓✓✓\n\n")
