#!/usr/bin/env Rscript

# =============================================================================
# 01_load_data_create_seurat.R
#
# Purpose:
#   - Recursively load raw 10X MTX matrices (12 samples) from GEO GSE231993
#   - Merge into one sparse count matrix (safe union of genes)
#   - Build Seurat object + attach metadata (condition/replicate/sample_id)
#   - Compute basic QC metrics (percent.mt, percent.ribo)
#   - Save Seurat object + clinical summary table
#
# Inputs (expected):
#   data/raw/GSE231993/**/matrix.mtx.gz
#   data/raw/GSE231993/**/barcodes.tsv.gz
#   data/raw/GSE231993/**/features.tsv.gz
#   optional: metadata/sample_sheet.csv (preferred GSM → condition mapping)
#
# Outputs:
#   results/01_quality_control/sessionInfo.txt
#   results/01_quality_control/tables/clinical_summary.csv
#   results/01_quality_control/objects/seurat_raw.rds
# =============================================================================

# =========================
# Libraries 
# =========================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(stringr)
})

# =========================
# Project root detection
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

raw_dir <- file.path(root_dir, "data", "raw")         
out_dir <- file.path(root_dir, "results", "01_quality_control")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# DATA SOURCE INFORMATION
# =========================
cat("\n=== DATASET INFORMATION ===\n")
cat("Dataset: Du et al. (2023) Nature Communications\n")
cat("GEO: GSE231993\n")
cat("Title: Selective oxidative protection in ulcerative colitis\n")
cat("Data type: RAW COUNTS (10X format - MTX files)\n")
cat("Samples: 12 (4 Ulcerative Colitis, 4 Heahtly Control, 4 Ulcerative Colitis-self-control)\n")
cat("Focus: Colonocytes in ulcerative colitis\n\n")

# =========================
# 1) LOAD DATA - 10X MTX
# =========================
cat("Loading raw count matrices...\n")

if (!dir.exists(raw_dir)) {
  stop("Expected raw data folder not found: ", raw_dir,
       "\nMake sure you extracted GEO tar into data/raw/ (e.g., data/raw/GSE231993/...)")
}

matrix_files <- list.files(
  path = raw_dir,
  pattern = "matrix\\.mtx\\.gz$",
  full.names = TRUE,
  recursive = TRUE # Ensures that it looks into all subfolders (can happen in mtx scRNA-seq files)
)

cat("Found", length(matrix_files), "matrix.mtx.gz files\n")
if (length(matrix_files) == 0) {
  stop("No matrix.mtx.gz files found under: ", raw_dir)
}

load_10x_sample <- function(mtx_file) {
  gsm_id <- str_extract(mtx_file, "GSM[0-9]+")
  if (is.na(gsm_id)) stop("Could not extract GSM ID from path: ", mtx_file)
  
  barcode_file <- gsub("matrix\\.mtx\\.gz$", "barcodes.tsv.gz", mtx_file)
  feature_file <- gsub("matrix\\.mtx\\.gz$", "features.tsv.gz", mtx_file)
  
  if (!file.exists(barcode_file)) stop("Missing barcodes file: ", barcode_file)
  if (!file.exists(feature_file)) stop("Missing features file: ", feature_file)
  
  cat("  Loading", gsm_id, "...\n")
  
  mat <- Seurat::ReadMtx(
    mtx = mtx_file,
    cells = barcode_file,
    features = feature_file,
    feature.column = 2
  )
  
  rownames(mat) <- make.unique(rownames(mat))
  colnames(mat) <- paste0(gsm_id, "_", colnames(mat))
  
  cat("    ", nrow(mat), "genes x", ncol(mat), "cells\n")
  mat
}

count_matrices <- list()
for (mtx_file in matrix_files) {
  gsm_id <- str_extract(mtx_file, "GSM[0-9]+")
  if (is.na(gsm_id)) stop("Could not extract GSM ID from path: ", mtx_file)
  count_matrices[[gsm_id]] <- load_10x_sample(mtx_file)
}

cat("\nMerging samples (safe union of genes)...\n")

merge_sparse_counts <- function(mats) {
  all_genes <- unique(unlist(lapply(mats, rownames))) 
  cat("Total unique genes across all samples:", length(all_genes), "\n")
  
  mats_aligned <- lapply(mats, function(m) {
    missing <- setdiff(all_genes, rownames(m))
    if (length(missing) > 0) {
      z <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        dims = c(length(missing), ncol(m)),
        dimnames = list(missing, colnames(m))
      )
      m <- rbind(m, z)
    }
    m[all_genes, , drop = FALSE]
  })
  
  Reduce(cbind, mats_aligned)
}

counts_merged <- merge_sparse_counts(count_matrices)

cat("✓ Data loaded\n")
cat("Genes:", nrow(counts_merged), "\n")
cat("Cells:", ncol(counts_merged), "\n")

rm(count_matrices) # Clearing memory
gc()

# =========================
# 2) PREPARE METADATA
# =========================
cat("\nPreparing metadata...\n")

cell_barcodes <- colnames(counts_merged)
metadata <- data.frame(
  Cell = cell_barcodes,
  stringsAsFactors = FALSE,
  row.names = cell_barcodes
)

metadata$gsm_id <- str_extract(metadata$Cell, "^GSM[0-9]+")

  
metadata <- metadata %>%
    mutate(
      condition = case_when(
        gsm_id %in% c("GSM7307094", "GSM7307095", "GSM7307096", "GSM7307097") ~ "UC Self-Control",
        gsm_id %in% c("GSM7307098", "GSM7307099", "GSM7307100", "GSM7307101") ~ "Ulcerative Colitis",
        gsm_id %in% c("GSM7307102", "GSM7307103", "GSM7307104", "GSM7307105") ~ "Healthy Control",
        TRUE ~ "Unknown"
      ),
      condition_short = case_when(
        condition == "Healthy Control" ~ "HC",
        condition == "UC Self-Control" ~ "UCSC",
        condition == "Ulcerative Colitis" ~ "UC",
        TRUE ~ "Unknown"
      ),
      replicate = case_when(
        gsm_id == "GSM7307094" ~ "rep1",
        gsm_id == "GSM7307095" ~ "rep2",
        gsm_id == "GSM7307096" ~ "rep3",
        gsm_id == "GSM7307097" ~ "rep4",
        gsm_id == "GSM7307098" ~ "rep1",
        gsm_id == "GSM7307099" ~ "rep2",
        gsm_id == "GSM7307100" ~ "rep3",
        gsm_id == "GSM7307101" ~ "rep4",
        gsm_id == "GSM7307102" ~ "rep5",
        gsm_id == "GSM7307103" ~ "rep6",
        gsm_id == "GSM7307104" ~ "rep7",
        gsm_id == "GSM7307105" ~ "rep8",
        TRUE ~ "unknown"
      ),
      sample_id = paste0(condition_short, "_", replicate)
    )

cat("✓ Metadata prepared\n")
cat("\nSample distribution:\n")
print(table(metadata$condition, useNA = "ifany"))
cat("\nCells per sample:\n")
print(table(metadata$sample_id, useNA = "ifany"))

# =========================
# 3) CREATE SEURAT OBJECT
# =========================
cat("\nCreating Seurat object...\n")

seurat_obj <- CreateSeuratObject(
  counts = counts_merged,
  project = "Ulcerative Colitis",
  min.cells = 3,
  min.features = 200,
  meta.data = metadata
)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

cat("✓ Seurat object created\n")
cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")

# Clinical / QC summary table
clinical_summary <- seurat_obj@meta.data %>%
  dplyr::group_by(condition, condition_short, replicate, sample_id) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_genes = round(mean(nFeature_RNA)),
    median_genes = median(nFeature_RNA),
    mean_umi = round(mean(nCount_RNA)),
    median_umi = median(nCount_RNA),
    mean_mito = round(mean(percent.mt), 2),
    .groups = "drop"
  )

print(clinical_summary)

write.csv(
  clinical_summary,
  file.path(table_dir, "clinical_summary.csv"),
  row.names = FALSE
)

# Save object for downstream scripts
saveRDS(seurat_obj, file.path(object_dir, "seurat_raw.rds"))

rm(counts_merged)
gc()

cat("\n✓✓✓ LOADING COMPLETE ✓✓✓\n")
cat("Saved: ", file.path(object_dir, "seurat_raw.rds"), "\n\n")
