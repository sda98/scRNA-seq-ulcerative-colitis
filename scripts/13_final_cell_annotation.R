#!/usr/bin/env Rscript

# =============================================================================
# 13_final_cell_annotation.R
#
# Purpose:
#   1) Assign final (manual / marker-based) cell type annotations to res=0.4 clusters
#   2) Run SingleR (Human Primary Cell Atlas) for automated validation
#   3) Remove the low-quality cluster
#   4) Add broad cell type categories (Immune/Epithelial/Stromal/Neural)
#   5) Save annotated Seurat objects
#
# Inputs (expected to already exist in the environment):
#   - seurat_integrated
#   - object_dir
#
# Outputs:
#   - seurat_integrated_with_annotations.rds
#   - seurat_clean_annotated.rds
#
# Notes:
#   - Uses Idents = "RNA_snn_res.0.4"
#   - Excludes: "Low-quality (excluded)"
# =============================================================================


# =============================================================================
# 17) CELL TYPE ANNOTATION
# =============================================================================

cat("\n=== CELL TYPE ANNOTATION ===\n")

# Set resolution 0.4 as active identity
Idents(seurat_integrated) <- "RNA_snn_res.0.4"

# Final annotations from marker-based analysis
final_annotations <- c(
  "0"  = "CD4+ T cells",
  "1"  = "Plasma cells (IgA+)",
  "2"  = "Cytotoxic lymphocytes (CD8+ T / NK)",
  "3"  = "B cells (naive/memory)",
  "4"  = "Low-quality (excluded)",
  "5"  = "Plasma cells (IgG+)",
  "6"  = "Colonocytes (less mature)",
  "7"  = "Colonocytes (more mature)",
  "8"  = "Inflammatory fibroblasts",
  "9"  = "Macrophages",
  "10" = "Cycling B cells",
  "11" = "Goblet cells",
  "12" = "Endothelial cells",
  "13" = "Mast cells",
  "14" = "Enteric glia",
  "15" = "Tuft cells",
  "16" = "Pericytes / SMCs",
  "17" = "Germinal center B cells",
  "18" = "Enteroendocrine cells (L-cell enriched)"
)

# Add annotations to metadata
cluster_ids  <- as.character(Idents(seurat_integrated))
final_labels <- final_annotations[cluster_ids]
names(final_labels) <- colnames(seurat_integrated)

seurat_integrated <- AddMetaData(
  seurat_integrated,
  metadata = final_labels,
  col.name = "celltype_final"
)

cat("✓ Final annotations added\n")
cat("\nCell type distribution:\n")
print(table(seurat_integrated$celltype_final, useNA = "ifany"))


# =============================================================================
# SingleR automated annotation (for validation)
# =============================================================================

cat("\n=== Running SingleR annotation ===\n")
cat("Reference: Human Primary Cell Atlas\n\n")

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

# Add SingleR results to Seurat object
seurat_integrated$singler_celltype        <- pred$labels
seurat_integrated$singler_celltype_pruned <- pred$pruned.labels
seurat_integrated$singler_score           <- apply(pred$scores, 1, max)

cat("✓ SingleR annotation complete\n")
cat("\nSingleR annotation distribution:\n")
print(table(seurat_integrated$singler_celltype))


# =============================================================================
# Step 2: Remove low-quality cluster
# =============================================================================

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


# =============================================================================
# Step 4: Create broad category annotations
# =============================================================================

cat("\nAdding broad cell type categories...\n")

broad_categories <- c(
  "CD4+ T cells"                           = "Immune - T cells",
  "Plasma cells (IgA+)"                    = "Immune - B/Plasma",
  "Cytotoxic lymphocytes (CD8+ T / NK)"    = "Immune - T cells",
  "B cells (naive/memory)"                 = "Immune - B/Plasma",
  "Plasma cells (IgG+)"                    = "Immune - B/Plasma",
  "Colonocytes*"                           = "Epithelial",
  "Colonocytes"                            = "Epithelial",
  "Inflammatory fibroblasts"               = "Stromal",
  "Macrophages"                            = "Immune - Myeloid",
  "Cycling B cells"                        = "Immune - B/Plasma",
  "Goblet cells"                           = "Epithelial",
  "Endothelial cells"                      = "Stromal",
  "Mast cells"                             = "Immune - Myeloid",
  "Enteric glia"                           = "Neural",
  "Tuft cells"                             = "Epithelial",
  "Pericytes / SMCs"                       = "Stromal",
  "Germinal center B cells"                = "Immune - B/Plasma",
  "Enteroendocrine cells (L-cell enriched)" = "Epithelial"
)

# Create named vector with cell barcodes
broad_labels <- broad_categories[seurat_clean$celltype_final]
names(broad_labels) <- colnames(seurat_clean)

seurat_clean <- AddMetaData(
  seurat_clean,
  metadata = broad_labels,
  col.name = "celltype_broad"
)

# Make it a factor with consistent order for plotting
seurat_clean$celltype_broad <- factor(
  seurat_clean$celltype_broad,
  levels = c(
    "Immune - T cells",
    "Immune - B/Plasma",
    "Immune - Myeloid",
    "Epithelial",
    "Stromal",
    "Neural"
  )
)

cat("✓ Broad categories added\n")
print(table(seurat_clean$celltype_broad))


# =============================================================================
# SAVE ANNOTATED OBJECTS (critical - do this once!)
# =============================================================================

cat("\n=== Saving annotated objects ===\n")

saveRDS(
  seurat_integrated,
  file.path(object_dir, "seurat_integrated_with_annotations.rds")
)
saveRDS(
  seurat_clean,
  file.path(object_dir, "seurat_clean_annotated.rds")
)

cat("✓ Annotated objects saved:\n")
cat("  - seurat_integrated_with_annotations.rds\n")
cat("  - seurat_clean_annotated.rds\n")


# =============================================================================
# Clean up
# =============================================================================

rm(seurat_obj, seurat_obj_filtered)
gc()
