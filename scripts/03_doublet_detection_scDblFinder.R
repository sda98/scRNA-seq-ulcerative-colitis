#!/usr/bin/env Rscript

# =============================================================================
# 03_doublet_detection_scDblFinder.R
#
# Purpose:
#   1) Load QC-filtered Seurat object from Script 02
#   2) Add readable per-sample labels (for plotting)
#   3) Run scDblFinder doublet detection per sample_id
#   4) Save doublet rates table + per-sample score histograms
#   5) Remove predicted doublets (keep singlets only)
#   6) Save singlet-only Seurat object for downstream scripts
#
# Inputs:
#   results/01_quality_control/objects/seurat_qc_filtered.rds
#
# Outputs:
#   results/01_quality_control/tables/doublet_rate_by_sample.csv
#   results/01_quality_control/plots/doublets_by_sample.png
#   results/01_quality_control/objects/seurat_singlets.rds
#
# Notes:
#   - Doublet calling is performed per sample (samples = "sample") to avoid
#     cross-sample effects.
# =============================================================================

# =========================
# Libraries (required)
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggh4x)
  library(SingleCellExperiment)
  library(scDblFinder)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "01_quality_control", "objects", "seurat_qc_filtered.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 02_qc_metrics_and_plots.R first.")
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
cat("\nLoading QC-filtered Seurat object...\n")
seurat_obj_filtered <- readRDS(in_obj)

# Basic sanity checks
if (!"sample_id" %in% colnames(seurat_obj_filtered@meta.data)) {
  stop("metadata column 'sample_id' not found. It should exist from Scripts 01/02.")
}
if (!"gsm_id" %in% colnames(seurat_obj_filtered@meta.data)) {
  stop("metadata column 'gsm_id' not found. It should exist from Script 01.")
}

# =============================================================================
# 1) ADD READABLE SAMPLE LABELS (FOR PLOTTING)
# =============================================================================
seurat_obj_filtered@meta.data <- seurat_obj_filtered@meta.data %>%
  mutate(
    sample_label = case_when(
      gsm_id == "GSM7307094" ~ "UC Self-Control 1",
      gsm_id == "GSM7307095" ~ "UC Self-Control 2",
      gsm_id == "GSM7307096" ~ "UC Self-Control 3",
      gsm_id == "GSM7307097" ~ "UC Self-Control 4",
      gsm_id == "GSM7307098" ~ "Ulcerative Colitis 1",
      gsm_id == "GSM7307099" ~ "Ulcerative Colitis 2",
      gsm_id == "GSM7307100" ~ "Ulcerative Colitis 3",
      gsm_id == "GSM7307101" ~ "Ulcerative Colitis 4",
      gsm_id == "GSM7307102" ~ "Healthy Control 5",
      gsm_id == "GSM7307103" ~ "Healthy Control 6",
      gsm_id == "GSM7307104" ~ "Healthy Control 7",
      gsm_id == "GSM7307105" ~ "Healthy Control 8",
      TRUE ~ sample_id
    )
  )

cat("✓ Sample labels added\n")

# =============================================================================
# 2) RUN scDblFinder (PER SAMPLE)
# =============================================================================
cat("\nRunning doublet detection (scDblFinder)...\n")

sce <- as.SingleCellExperiment(seurat_obj_filtered)
colData(sce)$sample <- seurat_obj_filtered$sample_id

sce <- scDblFinder(sce, samples = "sample")

seurat_obj_filtered$scDblFinder.class <- colData(sce)$scDblFinder.class
seurat_obj_filtered$scDblFinder.score <- colData(sce)$scDblFinder.score

cat("✓ scDblFinder complete\n")

# =============================================================================
# 3) EXPORT DOUBLET RATES (PER SAMPLE)
# =============================================================================
tab <- table(seurat_obj_filtered$sample_id, seurat_obj_filtered$scDblFinder.class)
dbl_rate_df <- as.data.frame.matrix(prop.table(tab, margin = 1))
dbl_rate_df$sample <- rownames(dbl_rate_df)

if (!"singlet" %in% colnames(dbl_rate_df)) dbl_rate_df$singlet <- NA_real_
if (!"doublet" %in% colnames(dbl_rate_df)) dbl_rate_df$doublet <- NA_real_

write.csv(
  dbl_rate_df[, c("sample", "singlet", "doublet")],
  file.path(table_dir, "doublet_rate_by_sample.csv"),
  row.names = FALSE
)

cat("\nDoublet summary (all samples):\n")
print(table(seurat_obj_filtered$scDblFinder.class))

# =============================================================================
# 4) PLOT DOUBLET SCORE DISTRIBUTIONS (FACET BY SAMPLE)
# =============================================================================
df_dbl_sample <- data.frame(
  score = colData(sce)$scDblFinder.score,
  class = factor(colData(sce)$scDblFinder.class, levels = c("singlet", "doublet")),
  sample_label = factor(seurat_obj_filtered$sample_label)
)

p_dbl_sample <- ggplot(df_dbl_sample, aes(x = score, fill = class)) +
  geom_histogram(
    bins = 50,
    color = "black",
    linewidth = 0.15,
    alpha = 0.85,
    position = "identity"
  ) +
  ggh4x::facet_wrap2(~ sample_label, scales = "fixed", axes = "all", ncol = 4) +
  labs(
    x = "scDblFinder score",
    y = "Number of cells",
    fill = "Class"
  ) +
  theme_classic(base_size = 32) +
  theme(
    strip.text = element_text(face = "bold", size = 27),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.title = element_text(face = "bold", size = 40),
    axis.title = element_text(face = "bold", size = 55),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0, unit = "pt")),
    axis.text = element_text(size = 28),
    legend.text = element_text(face = "bold", size = 35)
  )

ggsave(
  file.path(plot_dir, "doublets_by_sample.png"),
  p_dbl_sample, width = 24, height = 14, dpi = 600
)

cat("✓ Doublet histogram saved\n")

# =============================================================================
# 5) REMOVE DOUBLET CELLS (KEEP SINGLET ONLY)
# =============================================================================
cat("\nRemoving predicted doublets...\n")
seurat_obj_filtered <- subset(seurat_obj_filtered, subset = scDblFinder.class == "singlet")
cat("Cells after doublet removal:", ncol(seurat_obj_filtered), "\n")

rm(sce)
gc()

# =============================================================================
# 6) SAVE SINGLET-ONLY OBJECT
# =============================================================================
out_obj <- file.path(object_dir, "seurat_singlets.rds")
saveRDS(seurat_obj_filtered, out_obj)

cat("\n✓ Saved singlet-only object: ", out_obj, "\n")
cat("\n✓✓✓ DOUBLET DETECTION COMPLETE ✓✓✓\n\n")
