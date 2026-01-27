# scRNA-seq Analysis of Human Colon Biopsies in Ulcerative Colitis (GSE231993 / PRJNA970371)

This repository contains an end-to-end **single-cell RNA-seq (scRNA-seq)** workflow in **R** using **Seurat + Harmony** for a public UC colon biopsy dataset (Du et al., Nat Commun 2023).  
The pipeline loads raw **10X MTX** matrices, performs **QC + filtering (robust per-sample bounds)**, detects/removes **doublets (scDblFinder)**, runs **Harmony integration**, clusters cells, identifies markers, performs **cell type annotation (manual + SingleR validation)**, validates colonocyte maturation using **CytoTRACE2**, performs **cell-type composition analysis** (including **scCODA-ready** inputs), and runs **pseudobulk DESeq2** for targeted questions (e.g., **Colonocytes (less differentiated)**).

> Note: In the current “single-script” version, many outputs are written to `results/01_quality_control/` for convenience.  
> The planned repo structure below separates outputs by analysis stage (QC, integration, annotation, composition, pseudobulk, etc.).

---

## Study design (from GEO / BioProject)

- Organism: *Homo sapiens*
- Experiment type: Expression profiling by high throughput sequencing (scRNA-seq)
- Groups (12 samples):
  - **UCSC (UC self-control):** GSM7307094–GSM7307097
  - **UC inflamed:** GSM7307098–GSM7307101
  - **HC healthy controls:** GSM7307102–GSM7307105
- BioProject: **PRJNA970371**

---

## Data sources

### GEO series
- Dataset: **GSE231993**
- GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993

### Samples (12 total)
- **UCSC (UC self-control):** GSM7307094, GSM7307095, GSM7307096, GSM7307097
- **UC inflamed:** GSM7307098, GSM7307099, GSM7307100, GSM7307101
- **HC healthy controls:** GSM7307102, GSM7307103, GSM7307104, GSM7307105

### Download (recommended)
From the GEO page, download the supplementary archive:
- `GSE231993_RAW.tar` (10X MTX files)

Unpack into the repo (example):
```bash
mkdir -p data/raw/GSE231993
tar -xf GSE231993_RAW.tar -C data/raw/GSE231993


The loader expects each sample folder to contain the standard 10X files:
matrix.mtx.gz, features.tsv.gz (or genes.tsv.gz), barcodes.tsv.gz.

Repository structure 
data/
  raw/                         # (NOT committed) unpacked GSM folders with 10X MTX files
metadata/
  sample_sheet.csv             # optional: curated mapping GSM → condition/replicate
scripts/
  00_run_all.R
  01_load_data_create_seurat.R
  02_qc_metrics_and_plots.R
  03_qc_filtering_robust_bounds.R
  04_doublet_detection_scDblFinder.R
  05_preprocess_pca_horns_parallel.R
  06_harmony_integration_umap.R
  07_clustering_markers.R
  08_cluster_qc_and_cleanup.R
  09_celltype_annotation_manual.R
  10_singler_validation.R
  11_colonocyte_subset_cytotrace2.R
  12_composition_tables_and_plots.R
  13_prepare_sccoda_inputs.R
  14_run_sccoda_reticulate_or_python.py
  15_pseudobulk_prepare.R
  16_pseudobulk_deseq2.R
  17_summary_figures_for_readme.R
results/
  01_quality_control/
    plots/ tables/ objects/
  02_integration_clustering/
    plots/ tables/ objects/
  03_annotation/
    plots/ tables/ objects/
  04_cytotrace_colonocytes/
    plots/ tables/ objects/
  05_composition/
    plots/ tables/
  06_sccoda/
    plots/ tables/
  07_pseudobulk_deseq2/
    Colonocytes_less_differentiated/
      plots/ tables/

Scripts list (planned)
Orchestration

scripts/00_run_all.R
Runs the full workflow end-to-end by sourcing the scripts below in order.

Data loading + QC

scripts/01_load_data_create_seurat.R
Loads all 10X MTX samples, merges them, and creates the initial Seurat object with sample metadata.

scripts/02_qc_metrics_and_plots.R
Computes QC metrics and produces QC plots and per-sample summaries.

scripts/03_qc_filtering_robust_bounds.R
Computes robust per-sample QC bounds and filters low-quality cells.

scripts/04_doublet_detection_scDblFinder.R
Runs scDblFinder, visualizes doublets, removes predicted doublets, and exports summaries.

Preprocessing + integration + clustering

scripts/05_preprocess_pca_horns_parallel.R
Performs LogNormalize/HVG selection/ScaleData/PCA and Horn’s parallel analysis to choose the number of PCs.

scripts/06_harmony_integration_umap.R
Runs Harmony integration (batch = sample_id) and generates UMAPs before/after integration.

scripts/07_clustering_markers.R
Performs clustering at selected resolutions and identifies marker genes per cluster.

Cleanup + annotation + validation

scripts/08_cluster_qc_and_cleanup.R
Cluster-level QC diagnostics and removal of low-quality clusters when justified.

scripts/09_celltype_annotation_manual.R
Manual cluster-to-celltype mapping, creation of broad classes (Immune/Epithelial/Stromal), and UMAPs by annotation.

scripts/10_singler_validation.R
SingleR validation (Human Primary Cell Atlas via celldex) + manual vs SingleR cross-tabulation.

Colonocyte maturation

scripts/11_colonocyte_subset_cytotrace2.R
Subsets colonocytes, runs CytoTRACE2, and renames into less vs more differentiated colonocytes.

Composition + scCODA

scripts/12_composition_tables_and_plots.R
Cell counts/proportions per sample and summaries by condition.

scripts/13_prepare_sccoda_inputs.R
Exports scCODA-ready sccoda_counts.csv and sccoda_metadata.csv.

scripts/14_run_sccoda_reticulate_or_python.py
Optional: runs scCODA either as a standalone Python script or via reticulate.

Pseudobulk DESeq2

scripts/15_pseudobulk_prepare.R
Builds pseudobulk counts by sample for a chosen cell type and applies minimum-cell filters.

scripts/16_pseudobulk_deseq2.R
Runs DESeq2 contrasts (e.g., UC vs HC, UC vs UCSC) and exports plots/tables.

README figures

scripts/17_summary_figures_for_readme.R
Optional: creates a small set of “portfolio-ready” summary panels for the README.

Pipeline overview
1) Load raw 10X MTX data (12 samples)

Recursively finds matrix.mtx.gz files

Reads matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz

Prefixes barcodes with GSM..._ to ensure uniqueness

Merges into a single Seurat object

Adds sample metadata:

gsm_id

sample_id

condition (HC / UC / UCSC)

replicate

2) QC metrics and visualization

Computes:

percent.mt (MT- genes)

percent.ribo (RPS/RPL genes)

Saves per-sample QC summaries and QC plots

3) QC filtering (robust per-sample bounds)

Computes robust, sample-specific bounds for QC metrics

Filters low-quality cells using these bounds

Exports bounds table and “before vs after” QC plots

4) Doublet detection (scDblFinder)

Runs scDblFinder

Visualizes doublets by sample

Removes predicted doublets

Exports doublet-rate tables

5) Preprocessing + PCA + Horn’s parallel analysis

LogNormalize, FindVariableFeatures, ScaleData, RunPCA

Horn’s parallel analysis to estimate the number of informative PCs

Saves elbow plot and PCA diagnostics

6) Batch effect visualization + Harmony integration

UMAP before integration by sample

Harmony integration using sample_id as batch

UMAP after integration by sample

7) Clustering + markers

Clusters at chosen resolutions (e.g., 0.3, 0.4)

Identifies marker genes per cluster

Exports:

full markers table

top markers per cluster (top 5/10/20)

marker summary table

Generates cluster UMAPs (with and without labels)

8) Cluster QC and cleanup

Cluster-level QC plots (nUMI/nGene/mt/ribo)

Removes low-quality cluster(s) when justified

9) Cell type annotation (manual) + broad categories

Manual mapping produces celltype_final and celltype_broad (Immune / Epithelial / Stromal).
UMAPs are generated colored by:

cluster

broad cell type

final cell type

final cell type split by condition

10) SingleR validation

Runs SingleR using Human Primary Cell Atlas reference (via celldex)

Exports cross-tabulation (manual vs SingleR labels)

11) Colonocyte maturation validation (CytoTRACE2)

Subsets colonocytes and validates maturation ordering

Renames into:

Colonocytes (less differentiated)

Colonocytes (more differentiated)

Saves CytoTRACE2 plots and colonocyte marker DotPlots

12) Cell-type composition analysis + scCODA inputs

Counts and proportions per sample and cell type

Summaries by condition (HC / UC / UCSC)

Exports scCODA-ready:

sccoda_counts.csv

sccoda_metadata.csv

13) scCODA compositional modeling (optional)

Runs scCODA via reticulate or a Python script

Saves posterior / effect plots

14) Pseudobulk DESeq2 (targeted)

Builds pseudobulk counts by sample for a chosen cell type

Applies minimum-cell filter per sample (e.g., min_cells = 50)

Runs DESeq2 contrasts (example):

UC vs HC

UC vs UCSC

Saves PCA/MA/volcano + annotated DE tables

Expected outputs (examples)

Outputs are written under results/ and separated by analysis stage (planned).

Quality control (results/01_quality_control/)
Plots

QC violin plots (before/after)

QC scatter by sample

Doublets by sample

Tables

QC bounds table by sample

Per-sample QC summary

Doublet rate table by sample

Objects

Intermediate Seurat objects (filtered / singlet-only)

Integration + clustering (results/02_integration_clustering/)
Plots

UMAP before Harmony by sample

UMAP after Harmony by sample

Cluster UMAPs (labeled + clean/no-label versions)

PCA/Harmony embedding diagnostics

Tables

Markers tables (full + top N per cluster)

Marker summary per cluster

Objects

Integrated Seurat object + clustered object

Annotation (results/03_annotation/)
Plots

UMAP by celltype_broad

UMAP by celltype_final

UMAP by celltype_final split by condition

Tables

Manual vs SingleR cross-tabulation

Objects

Annotated Seurat object

CytoTRACE2 / colonocytes (results/04_cytotrace_colonocytes/)
Plots

DotPlot of colonocyte differentiation markers

CytoTRACE2 score plots confirming maturation ordering

Tables

Colonocyte subset summaries (if exported)

Composition (results/05_composition/)
Tables

Cell counts per sample

Cell type proportions per sample

Condition-level summaries

Fold-change summaries

scCODA (results/06_sccoda/)
Tables

sccoda_counts.csv

sccoda_metadata.csv

Plots

scCODA effect/posterior plots (if run)

Pseudobulk DESeq2 (results/07_pseudobulk_deseq2/)

Example subfolder:

Colonocytes_less_differentiated/

Plots

Dispersion plot

MA plots (unshrunk + shrunk)

PCA plots (PC1–PC2; PC3–PC4)

Volcano plots

Tables

Annotated DESeq2 result tables for contrasts (UC vs HC; UC vs UCSC)

Requirements
R

Key packages:

Seurat, harmony, tidyverse, ggplot2

scDblFinder, SingleCellExperiment, Matrix

SingleR, celldex

CytoTRACE2

DESeq2, ashr

reticulate (optional; for scCODA)

Python (optional; only for scCODA)

numpy, pandas, anndata, sccoda

The pipeline exports your session info (recommended) to a text file under results/.

How to run (planned script-based workflow)

Clone the repository and set the working directory to the repo root.

Place unpacked GEO 10X folders under:

data/raw/GSE231993/

Run the pipeline scripts in order:

source("scripts/01_load_data_create_seurat.R")
source("scripts/02_qc_metrics_and_plots.R")
source("scripts/03_qc_filtering_robust_bounds.R")
source("scripts/04_doublet_detection_scDblFinder.R")
source("scripts/05_preprocess_pca_horns_parallel.R")
source("scripts/06_harmony_integration_umap.R")
source("scripts/07_clustering_markers.R")
source("scripts/08_cluster_qc_and_cleanup.R")
source("scripts/09_celltype_annotation_manual.R")
source("scripts/10_singler_validation.R")
source("scripts/11_colonocyte_subset_cytotrace2.R")
source("scripts/12_composition_tables_and_plots.R")
source("scripts/13_prepare_sccoda_inputs.R")
# optional scCODA:
# run scripts/14_run_sccoda_reticulate_or_python.py (Python) OR a reticulate-based R wrapper
source("scripts/15_pseudobulk_prepare.R")
source("scripts/16_pseudobulk_deseq2.R")
source("scripts/17_summary_figures_for_readme.R")


Or run everything:

source("scripts/00_run_all.R")

Citation

Du J, Zhang J, Wang L, Wang X, et al.
Selective oxidative protection leads to tissue topological changes orchestrated by macrophage during ulcerative colitis.
Nature Communications (2023).

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993

BioProject: PRJNA970371

::contentReference[oaicite:0]{index=0}
