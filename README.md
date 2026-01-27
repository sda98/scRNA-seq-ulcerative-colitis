# scRNA-seq Analysis of Human Colon Biopsies in Ulcerative Colitis

**Dataset:** GSE231993 / PRJNA970371  
**Publication:** Du et al., *Nature Communications* (2023)

End-to-end single-cell RNA-seq workflow in R using **Seurat + Harmony** for a public UC colon biopsy dataset. The pipeline performs QC, batch correction, cell type annotation, compositional analysis, and pseudobulk differential expression.

---

## Overview

This pipeline:

- Loads raw **10X MTX** matrices (12 samples)
- Performs **QC + filtering** with robust per-sample bounds
- Detects and removes **doublets** (scDblFinder)
- Runs **Harmony integration** for batch correction
- Clusters cells and identifies **marker genes**
- Performs **cell type annotation** (manual + SingleR validation)
- Validates colonocyte maturation using **CytoTRACE2**
- Runs **compositional analysis** (scCODA-ready outputs)
- Performs **pseudobulk DESeq2** for differential expression

---

## Study Design

| Group | Description | Samples |
|-------|-------------|---------|
| **HC** | Healthy controls | GSM7307102–GSM7307105 |
| **UCSC** | UC self-control (non-inflamed) | GSM7307094–GSM7307097 |
| **UC** | UC inflamed | GSM7307098–GSM7307101 |

- **Organism:** *Homo sapiens*
- **Tissue:** Colon biopsies
- **Technology:** 10X Genomics scRNA-seq

---

## Repository Structure

```
├── data/
│   └── raw/                          # 10X MTX files (not committed)
├── metadata/
│   └── sample_sheet.csv              # GSM → condition mapping
├── scripts/
│   ├── 00_run_all.R                  # Orchestration script
│   ├── 01_load_data_create_seurat.R
│   ├── 02_qc_metrics_and_plots.R
│   ├── 03_qc_filtering_robust_bounds.R
│   ├── 04_doublet_detection_scDblFinder.R
│   ├── 05_preprocess_pca_horns_parallel.R
│   ├── 06_harmony_integration_umap.R
│   ├── 07_clustering_markers.R
│   ├── 08_cluster_qc_and_cleanup.R
│   ├── 09_celltype_annotation_manual.R
│   ├── 10_singler_validation.R
│   ├── 11_colonocyte_subset_cytotrace2.R
│   ├── 12_composition_tables_and_plots.R
│   ├── 13_prepare_sccoda_inputs.R
│   ├── 14_run_sccoda_reticulate_or_python.py
│   ├── 15_pseudobulk_prepare.R
│   ├── 16_pseudobulk_deseq2.R
│   └── 17_summary_figures_for_readme.R
└── results/
    ├── 01_quality_control/
    ├── 02_integration_clustering/
    ├── 03_annotation/
    ├── 04_cytotrace_colonocytes/
    ├── 05_composition/
    ├── 06_sccoda/
    └── 07_pseudobulk_deseq2/
```

---

## Pipeline Steps

### 1. Data Loading
- Recursively finds and loads `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`
- Prefixes barcodes with sample ID for uniqueness
- Merges into single Seurat object with metadata (condition, replicate)

### 2. Quality Control
- Computes `percent.mt` and `percent.ribo` metrics
- Applies **robust per-sample bounds** using median ± MAD
- Generates QC violin plots and scatter plots

### 3. Doublet Detection
- Runs **scDblFinder** on each sample
- Visualizes and removes predicted doublets
- Exports doublet rate summaries

### 4. Normalization & Dimensionality Reduction
- LogNormalize → FindVariableFeatures → ScaleData → PCA
- **Horn's parallel analysis** to select optimal number of PCs

### 5. Batch Correction & Integration
- **Harmony integration** using `sample_id` as batch variable
- UMAP visualization before/after correction

### 6. Clustering & Marker Identification
- Louvain clustering at multiple resolutions
- `FindAllMarkers` for cluster-specific genes
- Exports full markers table and top N per cluster

### 7. Cell Type Annotation
- **Manual annotation** based on canonical markers
- **SingleR validation** against Human Primary Cell Atlas
- Assigns `celltype_final` and `celltype_broad` (Immune/Epithelial/Stromal)

### 8. Colonocyte Maturation Analysis
- Subsets colonocyte populations
- **CytoTRACE2** to validate differentiation ordering
- Splits into "less differentiated" vs "more differentiated"

### 9. Compositional Analysis
- Cell counts and proportions per sample/condition
- Exports **scCODA-ready** CSV files for compositional modeling

### 10. Pseudobulk Differential Expression
- Aggregates counts by sample for target cell type
- **DESeq2** contrasts: UC vs HC, UC vs UCSC
- Generates MA plots, volcano plots, PCA

---

## Key Outputs

| Directory | Contents |
|-----------|----------|
| `results/01_quality_control/` | QC plots, filtering bounds, doublet summaries |
| `results/02_integration_clustering/` | UMAP plots, marker tables, Seurat objects |
| `results/03_annotation/` | Cell type UMAPs, SingleR cross-tabulation |
| `results/04_cytotrace_colonocytes/` | CytoTRACE2 scores, differentiation markers |
| `results/05_composition/` | Cell counts, proportions, fold-changes |
| `results/06_sccoda/` | scCODA input files, effect plots |
| `results/07_pseudobulk_deseq2/` | DE tables, volcano plots, PCA |

---

## Requirements

### R packages
```r
# Core
Seurat, harmony, tidyverse, ggplot2

# QC & preprocessing  
scDblFinder, SingleCellExperiment, Matrix

# Annotation
SingleR, celldex, CytoTRACE2

# Differential expression
DESeq2, ashr

# Optional
reticulate  # for scCODA
```

### Python (optional, for scCODA)
```
numpy, pandas, anndata, sccoda
```

---

## Usage

### Download data
```bash
# From GEO supplementary files
mkdir -p data/raw/GSE231993
tar -xf GSE231993_RAW.tar -C data/raw/GSE231993
```

### Run pipeline
```r
# Option 1: Run all steps
source("scripts/00_run_all.R")

# Option 2: Run individual scripts
source("scripts/01_load_data_create_seurat.R")
source("scripts/02_qc_metrics_and_plots.R")
# ... continue through pipeline
```

---

## Citation

> Du J, Zhang J, Wang L, Wang X, et al. Selective oxidative protection leads to tissue topological changes orchestrated by macrophage during ulcerative colitis. *Nature Communications* (2023).

**Data availability:**
- GEO: [GSE231993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993)
- BioProject: [PRJNA970371](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA970371)

---

## Author

Sergey Dadoyan, M.Sc.  
[LinkedIn](https://www.linkedin.com/in/sergey-dadoyan-msc/) | [GitHub](https://github.com/sda98)
