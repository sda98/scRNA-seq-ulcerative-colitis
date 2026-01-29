# !/usr/bin/env Rscript

# =============================================================================
# 10_singler_cell_annotation_validation.R
#
# Purpose:
#   1) Load Seurat object that already contains manual cell type labels (celltype_final).
#   2) Run SingleR using the Human Primary Cell Atlas reference (celldex) for validation.
#   3) Export:
#        - SingleR label distribution table
#        - Manual vs SingleR cross-tabulation
#        - Dominant SingleR label per Manual label (simple mapping table)
#   4) Save Seurat object with SingleR metadata added.
#
# Input:
#   results/02_clustering_analysis/objects/seurat_integrated_annotated.rds
#
# Outputs:
#   results/02_clustering_analysis/sessionInfo.txt
#   results/02_clustering_analysis/tables/singler_label_counts.csv
#   results/02_clustering_analysis/tables/manual_vs_singler_crosstab.csv
#   results/02_clustering_analysis/tables/dominant_singler_cluster_per_manual.csv
#   results/02_clustering_analysis/objects/seurat_integrated_annotated_singler.rds
#
# Notes:
#   - Uses RNA assay, log-normalized expression values from the "data" layer.
#   - Uses SerialParam() to run SingleR in a single process (most reproducible, and the safest option).
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "02_clustering_analysis", "objects",
                    "seurat_integrated_annotated.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 09_celltype_annotation_manual.R first (or ensure annotated object exists).")
}

out_dir    <- file.path(root_dir, "results", "02_clustering_analysis")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# 1) LOAD SEURAT OBJECT
# =========================
cat("\nLoading annotated Seurat object...\n")
seurat_integrated <- readRDS(in_obj)

if (!"celltype_final" %in% colnames(seurat_integrated@meta.data)) {
  stop("Missing 'celltype_final' in metadata. This script expects manual annotations already added.")
}

cat("✓ Object loaded\n")

# =============================================================================
# 2) SingleR with HUMAN ATLAS
# =============================================================================
cat("\n=== Running SingleR annotation ===\n")
cat("Reference: Human Primary Cell Atlas (celldex)\n\n")

DefaultAssay(seurat_integrated) <- "RNA"

# Get log-normalized data for SingleR
test_data <- GetAssayData(seurat_integrated, assay = "RNA", layer = "data")

# Load reference
ref <- celldex::HumanPrimaryCellAtlasData()

cat("Running SingleR (this may take a few minutes)...\n")

pred <- SingleR(
  test      = test_data,
  ref       = ref,
  labels    = ref$label.main,
  de.method = "wilcox",
  fine.tune = TRUE,
  BPPARAM   = BiocParallel::SerialParam()
)

# Add SingleR results to Seurat metadata
seurat_integrated$singler_celltype        <- pred$labels
seurat_integrated$singler_celltype_pruned <- pred$pruned.labels
seurat_integrated$singler_score           <- apply(pred$scores, 1, max)

cat("✓ SingleR annotation complete\n")

# Export SingleR label distribution
singler_counts <- as.data.frame(table(seurat_integrated$singler_celltype, useNA = "ifany"))
colnames(singler_counts) <- c("singler_celltype", "n_cells")
write.csv(singler_counts, file.path(table_dir, "singler_label_counts.csv"), row.names = FALSE)

cat("\nSingleR annotation distribution:\n")
print(table(seurat_integrated$singler_celltype, useNA = "ifany"))

# =============================================================================
# 3) MANUAL VS SINGLER COMPARISON
# =============================================================================
cat("\n=== Manual vs SingleR comparison ===\n")

cross_tab <- table(
  Manual  = seurat_integrated$celltype_final,
  SingleR = seurat_integrated$singler_celltype
)

write.csv(
  as.data.frame.matrix(cross_tab),
  file.path(table_dir, "manual_vs_singler_crosstab.csv"),
  row.names = TRUE
)

cat("\n=== Dominant SingleR label per Manual label ===\n")

dominant_map <- data.frame(
  Manual  = rownames(cross_tab),
  SingleR = apply(cross_tab, 1, function(x) names(which.max(x))),
  stringsAsFactors = FALSE
)

print(dominant_map, row.names = FALSE)

write.csv(
  dominant_map,
  file.path(table_dir, "dominant_singler_cluster_per_manual.csv"),
  row.names = FALSE
)

cat("✓ Cross-tabulation saved\n")

# =============================================================================
# 4) SAVE THE SEURAT OBJECT
# =============================================================================
out_obj <- file.path(object_dir, "seurat_integrated_annotated_singler.rds")
saveRDS(seurat_integrated, out_obj)

cat("\n✓ Saved object with SingleR metadata: ", out_obj, "\n")
cat("\n✓✓✓ SingleR validation complete ✓✓✓\n\n")
