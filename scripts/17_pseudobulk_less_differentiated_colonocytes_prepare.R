#!/usr/bin/env Rscript

# =============================================================================
# 17_pseudobulk_less_mature_colonocytes_prepare.R
#
# Purpose:
#   1) Load cleaned/annotated Seurat object (Script 13 output).
#   2) Subset a target cell type (default: Colonocytes (less differentiated)).
#   3) Build pseudobulk RNA count matrix (genes x samples) by summing counts per sample.
#   4) Create DESeq2 dataset (dds) with condition design and basic gene filtering.
#   5) Run DESeq2 and save core QC output (dispersion plot) + intermediate objects.
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_clean.rds
#     - meta.data must include:
#         - celltype_final
#         - sample_id (or gsm_id)
#         - condition (or group)
#     - RNA assay must include: layer='counts'
#
# Outputs (under results/04_pseudobulk_deseq2):
#   sessionInfo.txt
#   objects/pseudobulk_<celltype>.rds            (list: pseudobulk matrix + sample_md)
#   objects/dds_<celltype>.rds                  (DESeq2 object after DESeq())
#   plots/pseudobulk_DESeq2/<celltype>/Dispersion_plot.png
#
# Notes:
#   - This script PREPARES + RUNS DESeq2 for the selected cell type.
#   - Downstream result extraction, contrasts, shrinkage, PCA, annotation can be a next script.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(tibble)
  library(DESeq2)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(
  root_dir, "results", "02_clustering_analysis", "objects",
  "seurat_integrated_clean.rds"
)
if (!file.exists(in_obj)) {
  stop(
    "Input object not found: ", in_obj, "\n",
    "Expected cleaned annotated object from Script 13."
  )
}

out_dir    <- file.path(root_dir, "results", "04_pseudobulk_deseq2")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(out_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  PSEUDOBULK DESEQ2 (PREP) — TARGET CELL TYPE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Input object :", in_obj, "\n")
cat("Output dir   :", out_dir, "\n\n")

# =========================
# Step 1: Load object
# =========================
cat("Loading Seurat object...\n")
seurat_clean <- readRDS(in_obj)
cat("✓ Object loaded\n")

# -------------------------
# Metadata column mapping
# -------------------------
md <- seurat_clean@meta.data

if (!"celltype_final" %in% colnames(md)) stop("Missing metadata column: celltype_final")

sample_col <- if ("sample_id" %in% colnames(md)) {
  "sample_id"
} else if ("gsm_id" %in% colnames(md)) {
  "gsm_id"
} else {
  stop("No sample identifier found (expected 'sample_id' or 'gsm_id').")
}

condition_col <- if ("condition" %in% colnames(md)) {
  "condition"
} else if ("group" %in% colnames(md)) {
  "group"
} else {
  stop("No condition column found (expected 'condition' or 'group').")
}

cat("Sample column   :", sample_col, "\n")
cat("Condition column:", condition_col, "\n\n")

# =========================
# Step 2: Settings
# =========================
ct <- "Colonocytes (less mature)" # Can choose any other cluster here 
min_cells <- 50

safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)

plot_dir_ct  <- file.path(plot_dir, "pseudobulk_DESeq2", safe_ct)
table_dir_ct <- file.path(table_dir, "pseudobulk_DESeq2", safe_ct)

dir.create(plot_dir_ct,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir_ct, recursive = TRUE, showWarnings = FALSE)

cat(paste(rep("-", 80), collapse = ""), "\n")
cat("Cell type :", ct, "\n")
cat("min_cells :", min_cells, "\n")
cat("Plot dir  :", plot_dir_ct, "\n")
cat("Table dir :", table_dir_ct, "\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# =========================
# Step 3: Create pseudobulk counts
# =========================
cat("Creating pseudobulk counts...\n")

cells_ct <- rownames(md)[md$celltype_final == ct]
cat("Cells in target cell type:", length(cells_ct), "\n")
if (length(cells_ct) == 0) stop("No cells found for celltype_final == '", ct, "'")

cells_per_sample <- md[cells_ct, , drop = FALSE] %>%
  dplyr::count(sample_id = .data[[sample_col]], name = "n_cells") %>%
  dplyr::arrange(dplyr::desc(n_cells))

print(cells_per_sample)

good_samples <- cells_per_sample$sample_id[cells_per_sample$n_cells >= min_cells]
cat("Samples with >=", min_cells, "cells:", length(good_samples), "\n")
if (length(good_samples) < 2) {
  stop("Not enough samples with >= min_cells to build pseudobulk (need at least 2).")
}

cells_use <- cells_ct[md[cells_ct, sample_col] %in% good_samples]
cat("Cells retained after sample filter:", length(cells_use), "\n\n")

# Raw counts (genes x cells)
counts_mat <- GetAssayData(seurat_clean, assay = "RNA", layer = "counts")[, cells_use]

# Aggregate by sample (genes x samples)
sample_ids <- factor(md[cells_use, sample_col], levels = good_samples)

pseudobulk <- sapply(levels(sample_ids), function(s) {
  Matrix::rowSums(counts_mat[, sample_ids == s, drop = FALSE])
})
pseudobulk <- as.matrix(pseudobulk)

cat("Pseudobulk matrix:", nrow(pseudobulk), "genes x", ncol(pseudobulk), "samples\n")

# =========================
# Step 4: Sample metadata aligned to pseudobulk
# =========================
cat("\nBuilding sample metadata...\n")

sample_md <- md %>%
  tibble::rownames_to_column(var = ".cell") %>%  # removes rownames safely
  dplyr::distinct(
    sample_id = .data[[sample_col]],
    condition = .data[[condition_col]]
  ) %>%
  dplyr::filter(sample_id %in% colnames(pseudobulk)) %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c("Healthy Control", "UC Self-Control", "Ulcerative Colitis")
    )
  ) %>%
  dplyr::arrange(match(sample_id, colnames(pseudobulk))) %>%
  tibble::column_to_rownames("sample_id")


# Verify alignment
stopifnot(identical(rownames(sample_md), colnames(pseudobulk)))
print(table(sample_md$condition))

# Save pseudobulk inputs as an RDS bundle
pseudobulk_bundle <- list(
  celltype     = ct,
  min_cells    = min_cells,
  pseudobulk   = pseudobulk,
  sample_md    = sample_md,
  sample_col   = sample_col,
  condition_col= condition_col
)
pb_path <- file.path(object_dir, paste0("pseudobulk_", safe_ct, ".rds"))
saveRDS(pseudobulk_bundle, pb_path)
cat("✓ Saved pseudobulk bundle:", pb_path, "\n")

# =========================
# Step 5: Create DESeq2 dataset (dds)
# =========================
cat("\nCreating DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(pseudobulk),
  colData   = sample_md,
  design    = ~ condition
)

# Filter lowly expressed genes (basic, conservative)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# =========================
# Step 6: Run DESeq2
# =========================
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)
cat("✓ DESeq2 finished\n")

dds_path <- file.path(object_dir, paste0("dds_", safe_ct, ".rds"))
saveRDS(dds, dds_path)
cat("✓ Saved dds object:", dds_path, "\n")

# =========================
# Step 7: QC plot (Dispersion)
# =========================
cat("\nSaving dispersion plot...\n")

png(
  filename = file.path(plot_dir_ct, "Dispersion_plot.png"),
  width = 8, height = 7, units = "in", res = 600
)
plotDispEsts(dds)
dev.off()

cat("✓ Dispersion plot saved\n")

cat("\n✓✓✓ PSEUDOBULK PREP COMPLETE ✓✓✓\n\n")

