# Single-Cell RNA-seq Analysis of Ulcerative Colitis Colon Biopsies

**Dataset:** GSE231993 / PRJNA970371  
**Publication:** Du et al., *Nature Communications* (2023)  (https://doi.org/10.1038/s41467-023-39173-2)

**Title:** "Selective oxidative protection leads to tissue topological changes orchestrated by macrophage during ulcerative colitis"

## Overview

This repository contains a comprehensive end-to-end single-cell RNA-seq analysis pipeline for studying colonocyte biology and immune cell populations in ulcerative colitis (UC). The workflow processes 12 colon biopsy samples using Seurat, Harmony, and DESeq2 to perform:

- Quality control with robust per-sample thresholding
- Doublet detection and removal
- Batch correction via Harmony integration
- Cell type annotation (manual + automated validation)
- Colonocyte differentiation state validation using CytoTRACE2
- Compositional analysis with scCODA
- Pseudobulk differential expression analysis

---

## Study Design

| **Group** | **Description** | **Samples** | **n** |
|-----------|----------------|-------------|-------|
| **HC** | Healthy Control | GSM7307102–GSM7307105 | 4 |
| **UCSC** | UC Self-Control (non-inflamed tissue from UC patients) | GSM7307094–GSM7307097 | 4 |
| **UC** | Ulcerative Colitis (inflamed) | GSM7307098–GSM7307101 | 4 |

- **Organism:** *Homo sapiens*
- **Tissue:** Colon biopsies
- **Technology:** 10X Genomics Chromium Single Cell 3' RNA-seq
- **Sequencing:** Illumina NovaSeq 6000

---

## Repository Structure
```
.
├── data/
│   └── raw/                           # 10X MTX files (user must download/extract)
├── scripts/
│   ├── 01_load_data_create_seurat.R
│   ├── 02_qc_metrics_and_plots.R
│   ├── 03_doublet_detection_scDblFinder.R
│   ├── 04_normalization_finding_most_variable_features.R
│   ├── 05_preprocess_pca_horns_parallel.R
│   ├── 06_harmony_integration_umap.R
│   ├── 07_clustering.R
│   ├── 08_clustering_finding_markers.R
│   ├── 09_celltype_annotation_manual.R
│   ├── 10_singler_cell_annotation_validation.R
│   ├── 11_clustering_qc_cluster4_investigation.R
│   ├── 12_colonocyte_subset_cytotrace2.R
│   ├── 13_final_cell_annotation.R
│   ├── 14_composition_tables_and_plots.R
│   ├── 15_prepare_sccoda_inputs.R
│   ├── 16_run_sccoda_reticulate.R
│   ├── 17_pseudobulk_less_mature_colonocytes_prepare.R
│   └── 18_pseudobulk_less_mature_colonocytes.R
└── results/
    ├── 01_quality_control/
    ├── 02_clustering_analysis/
    ├── 03_compositional_analysis/
    └── 04_pseudobulk_deseq2/
```

---

## Installation & Requirements

### R Packages
```r
# Core workflow
install.packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork"))

# Integration & batch correction
install.packages("harmony")

# Quality control
BiocManager::install(c("scDblFinder", "SingleCellExperiment"))

# Cell type annotation
BiocManager::install(c("SingleR", "celldex"))

# Differentiation analysis
install.packages("devtools")
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")

# Differential expression
BiocManager::install(c("DESeq2", "ashr", "org.Hs.eg.db"))

# Visualization & utilities
install.packages(c("ggrepel", "ggh4x", "ggforce", "viridis", "scales"))

# For scCODA integration
install.packages("reticulate")
```

### Python Environment (for scCODA)
```bash
# Create conda environment
conda create -n sccoda python=3.9
conda activate sccoda

# Install scCODA
pip install sccoda

# Additional dependencies
pip install numpy pandas anndata matplotlib
```

---

## Data Download & Preparation

### **CRITICAL STEP:** Download Raw Data

Before running any scripts, you **must** download and extract the raw data from GEO:

1. **Download the raw data archive:**
   - Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993
   - Download: `GSE231993_RAW.tar` (under "Supplementary file")

2. **Extract to the correct location:**
```bash
   # Create data directory
   mkdir -p data/raw
   
   # Extract the tar archive
   tar -xf GSE231993_RAW.tar -C data/raw/
```

3. **Verify the structure:**
```bash
   # You should see 12 sample subdirectories, each containing:
   # - matrix.mtx.gz
   # - barcodes.tsv.gz
   # - features.tsv.gz
   
   ls data/raw/
   # Expected output:
   # GSM7307094  GSM7307095  GSM7307096  ... GSM7307105
```

**Note:** The pipeline recursively searches `data/raw/` for all `matrix.mtx.gz` files, so the exact subdirectory structure doesn't matter as long as all files are under `data/raw/`.

---

## Pipeline Workflow

### **Script 01:** Data Loading
- Recursively finds all 10X MTX files under `data/raw/`
- Loads matrix, barcodes, and features for each sample
- Merges samples 
- Creates initial Seurat object with metadata (condition, replicate, sample_id)
- Computes basic QC metrics (percent.mt (% mitochondrial reads), percent.ribo (% ribosomal reads))
- **Output:** `seurat_raw.rds`, clinical summary table

### **Script 02:** Quality Control & Filtering
- Generates QC violin plots (nCount_RNA, nFeature_RNA, percent.mt)
- Computes **robust per-sample QC bounds** using MAD (Median Absolute Deviation)
- Applies per-sample lower bounds for nCount/nFeature + global mitochondrial cutoff (15%)
- Visualizes cells before/after filtering (bar plots, scatter plots)
- **Output:** `seurat_qc_filtered.rds`, QC plots, per-sample threshold table

### **Script 03:** Doublet Detection
- Runs scDblFinder on each sample independently (avoids cross-sample artifacts)
- Generates per-sample doublet score histograms
- Removes predicted doublets (keeps singlets only)
- **Output:** `seurat_singlets.rds`, doublet rate table

### **Script 04:** Normalization & Variable Features
- LogNormalize with scale factor 10,000
- Identifies 2,000 highly variable genes using `vst` method
- Generates styled scatter plot of mean expression vs. standardized variance
- **Output:** `seurat_norm_hvg.rds`, variable features plot

### **Script 05:** PCA & Dimensionality Selection
- Scales data using HVGs
- Runs PCA on scaled data
- **Horn's Parallel Analysis** to estimate optimal number of PCs (avoids arbitrary cutoffs)
- Generates elbow plot with recommended PC cutoff
- **Output:** `seurat_pca.rds`, elbow plot with PC recommendation

### **Script 06:** Harmony Integration & UMAP
- Visualizes batch effects **before integration** (UMAP colored by sample)
- Runs Harmony batch correction using `sample_id` as batch variable
- Performs standard workflow on Harmony embeddings:
  - UMAP projection
  - Nearest neighbor graph
  - Clustering at multiple resolutions (0.05–1.0)
- Sets default resolution to 0.4
- Visualizes batch effects **after integration** (UMAP colored by sample)
- **Output:** `seurat_harmony_umap.rds`, before/after integration UMAPs

### **Script 07:** Clustering Visualization
- Generates publication-style UMAP plots:
  - Resolution 0.4 with clean repelled cluster labels
  - Resolution 0.3 with dashed circle highlighting colonocyte region
  - Resolution 0.4 with dashed circle highlighting colonocyte region
- **Output:** Cluster UMAPs, `seurat_integrated_clustered.rds`

### **Script 08:** Marker Gene Identification
- Runs FindAllMarkers (Wilcoxon test) for all clusters at resolution 0.4
- Filters out mitochondrial, ribosomal, and hemoglobin genes
- Exports marker tables:
  - All significant markers
  - Top 5, 10, 20 markers per cluster
  - Per-cluster summary statistics
- **Output:** Marker gene tables (CSV)

### **Script 09:** Manual Cell Type Annotation
- Assigns biologically meaningful labels to each cluster based on canonical markers
- Creates DotPlot showing colonocyte differentiation markers grouped by category:
  - Pan-epithelial markers (EPCAM, LGALS4, etc.)
  - Less differentiated epithelial markers (SLC12A2, EPHB2, OLFM4, MKI67, etc.)
  - More differentiated epithelial markers (GUCA2A/B, SLC26A3, CA4, BEST4, KRT20)
- Exports cell type distribution table
- **Output:** `seurat_integrated_annotated.rds`, cell type counts, differentiation marker DotPlot

### **Script 10:** SingleR Machine Learning Validation
- Validates manual annotations using SingleR machine learning-based tool with Human Primary Cell Atlas reference
- Generates cross-tabulation of manual vs. SingleR labels
- Identifies dominant SingleR label for each manual cell type
- **Output:** `seurat_integrated_annotated_singler.rds`, validation tables

### **Script 11:** Cluster 4 QC Investigation
- Investigates whether Cluster 4 represents low-quality cells
- Performs detailed QC analysis:
  - Per-cluster QC metric summaries
  - Top 20 marker categorization (immune, cytoskeletal, translational, metabolic)
  - Sample composition analysis (pie chart)
  - **Sample-level paired comparison** (Cluster 4 vs. Others)
    - Uses per-sample medians (defensible approach, avoids pseudoreplication)
    - Paired Wilcoxon test on sample-level data
- **Output:** QC tables, marker plots, composition pie chart, bar comparison with p-values

### **Script 12:** CytoTRACE2 Differentiation Validation
- Subsets colonocyte populations (Clusters 6 & 7)
- Runs CytoTRACE2 to obtain independent differentiation scores
- Validates that Cluster 6 is less differentiated than Cluster 7:
  - Cell-level summary statistics
  - **Sample-level paired analysis** (per-sample medians)
  - Paired Wilcoxon test across samples
- Generates visualizations:
  - Boxplot with significance annotation
  - Spaghetti plot showing per-sample consistency
- **Output:** `seurat_colonocytes_cytotrace2.rds`, CytoTRACE2 tables & plots

### **Script 13:** Final Cell Type Annotation
- Removes low-quality cluster (Cluster 4)
- Adds broad cell type categories:
  - Immune - T cells
  - Immune - B/Plasma
  - Immune - Myeloid
  - Epithelial
  - Stromal
  - Neural
- Generates publication-ready UMAPs:
  - Final cell type labels with color-coded categories
  - Broad categories
  - Cell types split by condition (HC, UCSC, UC) with colonocyte highlight
- **Output:** `seurat_integrated_clean.rds`, final annotation UMAPs

### **Script 14:** Compositional Summary
- Calculates cell counts and proportions per sample
- Generates per-condition summary statistics (mean, SD, median)
- **Output:** Cell count matrix, composition summary table

### **Script 15:** scCODA Input Preparation
- Creates scCODA-compatible input files:
  - `sccoda_counts.csv` (samples × cell types count matrix)
  - `sccoda_metadata.csv` (sample_id, condition, condition_code)
- **Output:** scCODA input CSVs

### **Script 16:** scCODA Bayesian Compositional Analysis
- Runs scCODA via Python/reticulate
- Performs Bayesian modeling with condition as covariate (HC baseline)
- Extracts inclusion probabilities and effect sizes
- Generates per-cell-type composition boxplots with scCODA significance stars:
  - `*` = inclusion probability ≥ 0.7
  - `**` = inclusion probability ≥ 0.9
- **Output:** scCODA model summary, per-cell-type box plots with statistical annotations 

### **Script 17:** Pseudobulk DESeq2 Preparation
- Subsets target cell type (default: "Colonocytes (less mature)")
- Aggregates raw counts by sample (pseudobulk approach)
- Creates DESeq2 dataset with condition design
- Applies basic gene filtering (≥10 counts in ≥3 samples)
- Runs DESeq2 normalization and dispersion estimation
- **Output:** `pseudobulk_<celltype>.rds`, `dds_<celltype>.rds`, dispersion plot

### **Script 18:** Pseudobulk Differential Expression Analysis
- Performs PCA on variance-stabilized counts (PC1/2 and PC3/4)
- Differential expression contrasts:
  - **UC vs. Healthy Control**
  - **UC vs. Self-Control**
- For each contrast:
  - Extracts results with lfcShrink (ashr method)
  - Annotates genes using org.Hs.eg.db
  - Generates MA plots (unshrunk & shrunk)
  - Creates volcano plots with functional category labels
- **Output:** Annotated DE tables, PCA plots, MA plots, volcano plots

---

## Key Outputs

### Quality Control (`results/01_quality_control/`)
- QC violin plots (before/after filtering)
- Scatter plots with per-sample threshold lines
- Doublet detection histograms
- Variable features plot
- Elbow plot with recommended PCs

### Clustering & Integration (`results/02_clustering_analysis/`)
- UMAP plots (before/after Harmony integration)
- Cluster UMAP plots (multiple resolutions)
- Marker gene tables (all markers, top 5/10/20)
- Cell type annotation UMAPs
- SingleR validation cross-tabulation
- Cluster 4 QC investigation plots
- CytoTRACE2 differentiation validation plots
- Final cell type UMAPs (by condition)

### Compositional Analysis (`results/03_compositional_analysis/`)
- Cell count matrix (samples × cell types)
- Composition summary by condition
- scCODA input files
- Per-cell-type composition boxplots with significance

### Differential Expression (`results/04_pseudobulk_deseq2/`)
- PCA plots (PC1/2, PC3/4)
- Annotated DE tables (UC vs HC, UC vs UCSC)
- MA plots (unshrunk & shrunk)
- Volcano plots with functional annotations
- Dispersion estimation plots

---

## Cell Type Annotations

| **Cluster** | **Cell Type** | **Category** |
|-------------|---------------|--------------|
| 0 | CD4+ T cells | Immune - T cells |
| 1 | Plasma cells (IgA+) | Immune - B/Plasma |
| 2 | Cytotoxic lymphocytes (CD8+ T / NK) | Immune - T cells |
| 3 | B cells (naive/memory) | Immune - B/Plasma |
| 4 | Low-quality (excluded) | — |
| 5 | Plasma cells (IgG+) | Immune - B/Plasma |
| 6 | Colonocytes (less mature) | Epithelial |
| 7 | Colonocytes (more mature) | Epithelial |
| 8 | Inflammatory fibroblasts | Stromal |
| 9 | Macrophages | Immune - Myeloid |
| 10 | Cycling B cells | Immune - B/Plasma |
| 11 | Goblet cells | Epithelial |
| 12 | Endothelial cells | Stromal |
| 13 | Mast cells | Immune - Myeloid |
| 14 | Enteric glia | Neural |
| 15 | Tuft cells | Epithelial |
| 16 | Pericytes / SMCs | Stromal |
| 17 | Germinal center B cells | Immune - B/Plasma |
| 18 | Enteroendocrine cells (L-cell enriched) | Epithelial |

---

## Usage

### 1. Download and Prepare Data
See **Data Download & Preparation** section above.

### 2. Downoload or Copy this Repository and Run Pipeline Sequentially
```r
# Set working directory to project root
setwd("/path/to/UC_scRNAseq_analysis")

# Run scripts in order
source("scripts/01_load_data_create_seurat.R")
source("scripts/02_qc_metrics_and_plots.R")
source("scripts/03_doublet_detection_scDblFinder.R")
source("scripts/04_normalization_finding_most_variable_features.R")
source("scripts/05_preprocess_pca_horns_parallel.R")
source("scripts/06_harmony_integration_umap.R")
source("scripts/07_clustering.R")
source("scripts/08_clustering_finding_markers.R")
source("scripts/09_celltype_annotation_manual.R")
source("scripts/10_singler_cell_annotation_validation.R")
source("scripts/11_clustering_qc_cluster4_investigation.R")
source("scripts/12_colonocyte_subset_cytotrace2.R")
source("scripts/13_final_cell_annotation.R")
source("scripts/14_composition_tables_and_plots.R")
source("scripts/15_prepare_sccoda_inputs.R")
source("scripts/16_run_sccoda_reticulate.R")
source("scripts/17_pseudobulk_less_mature_colonocytes_prepare.R")
source("scripts/18_pseudobulk_less_mature_colonocytes.R")
```

---

## Citation

**Original Publication:**  
Du J, Zhang J, Wang L, Wang X, et al. (2023). Selective oxidative protection leads to tissue topological changes orchestrated by macrophage during ulcerative colitis. *Nature Communications*, 14, 4868. https://doi.org/10.1038/s41467-023-40568-w

**Data Availability:**
- GEO: [GSE231993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993)
- BioProject: [PRJNA970371](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA970371)

---

## Requirements

- **R** (tested with **R 4.5.2**)  
- Tested on **Windows 11** using **RStudio 2026.01.0 Build 392**

---

## Author

**Sergey Dadoyan, M.Sc.**  
Bioinformatician 
[LinkedIn](https://linkedin.com/in/sergeydadoyan) | [GitHub](https://github.com/sergeydadoyan)

---

## License

This analysis code is provided for research and educational purposes. Please cite both the original publication and this repository when using this dataset or analysis approach.
