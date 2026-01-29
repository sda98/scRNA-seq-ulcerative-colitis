##!/usr/bin/env Rscript

# =============================================================================
# 16_run_sccoda_reticulate.R
#
# Purpose:
#   Run scCODA Bayesian compositional analysis (via Python) using R + reticulate:
#     1) Configure/verify the Python environment from R (reticulate).
#     2) Load scCODA inputs (counts + metadata) and build an AnnData object.
#     3) Fit scCODA model (HMC sampling) with condition_code (HC baseline).
#     4) Extract scCODA results (effect_df, summary).
#     5) Create per-celltype composition boxplots with scCODA inclusion-probability
#        star brackets, and save one PNG per cell type.
#
# Inputs (created by Script 15_prepare_sccoda_inputs.R):
#   - sccoda_counts.csv    (samples x celltypes; counts)
#   - sccoda_metadata.csv  (sample_id + condition_code, etc.)
#   NOTES 
#   - This script currently uses hard-coded Windows paths inside py_run_string().
#        
#
# Outputs:
#   - scCODA model summary printed to console
#   - A folder of PNGs (one per cell type), saved to the output directory defined
#     in the Python plotting block (out_dir_main; fallback out_dir_fb).
#
# Requirements / Setup (IMPORTANT):
#   - This is an R script that *calls Python* through reticulate.
#   - You must create/activate a working Python environment BEFORE running:
#       * Python >= 3.8
#       * scCODA installed: pip install sccoda
#       * scCODA depends on tensorflow (>=2.4) and tensorflow-probability (>=0.12);
#         GPU versions are not recommended/tested for scCODA. :contentReference[oaicite:0]{index=0}
#     Also required by your code: numpy, pandas, anndata, matplotlib.
#
# Notes:
#   - In reticulate, set the Python env *before* calling py_config() / imports
#     (e.g., reticulate::use_condaenv("sccoda", required = TRUE)).
#   - condition_code is treated as categorical with levels: HC, UCSC, UC (HC baseline).
#   - Stars indicate scCODA inclusion probability thresholds in your plotting code:
#       *  IP >= 0.7,  ** IP >= 0.9
#   -  The current version of this script has been produced by S. Dadoyan 
#      with the assistance of Claude (Opus 4.5) and GPT-5.2 LLM models. The code was fully tested and customized independently by S. Dadoyan
#
# Reference:
#   Büttner et al. (2021) Nature Communications; scCODA (theislab).
# =============================================================================
# =============================================================================
# scCODA ANALYSIS (R + reticulate; Python backend)
# =============================================================================
# =========================
# Libraries
# =========================

suppressPackageStartupMessages({
  library(reticulate)
})

# Paths and I/O (scCODA)

root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

comp_dir  <- file.path(root_dir, "results", "03_compositional_analysis")
table_dir <- file.path(comp_dir, "tables")
plot_dir  <- file.path(comp_dir, "plots", "sccoda_celltype_plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

counts_csv <- file.path(table_dir, "sccoda_counts.csv")
meta_csv   <- file.path(table_dir, "sccoda_metadata.csv")

if (!file.exists(counts_csv)) {
  stop("Missing scCODA counts file: ", counts_csv, "\nRun 15_prepare_sccoda_inputs.R first.")
}
if (!file.exists(meta_csv)) {
  stop("Missing scCODA metadata file: ", meta_csv, "\nRun 15_prepare_sccoda_inputs.R first.")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  scCODA INPUTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Counts CSV :", counts_csv, "\n")
cat("Meta CSV   :", meta_csv, "\n")
cat("Plot dir   :", plot_dir, "\n\n")

# Make paths Python-friendly
counts_csv_py <- normalizePath(counts_csv, winslash = "/", mustWork = TRUE)
meta_csv_py   <- normalizePath(meta_csv,   winslash = "/", mustWork = TRUE)
plot_dir_py   <- normalizePath(plot_dir,   winslash = "/", mustWork = TRUE)

# =============================================================================
# Create AnnData object
# =============================================================================
cat("Creating AnnData object...\n")

py_run_string(sprintf('
counts_path = r"%s"
meta_path   = r"%s"
out_dir_str = r"%s"
', counts_csv_py, meta_csv_py, plot_dir_py))

py_run_string('
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path

out_dir = Path(out_dir_str)
out_dir.mkdir(parents=True, exist_ok=True)

print("Counts path:", counts_path)
print("Meta path  :", meta_path)
print("Out dir    :", str(out_dir))

counts_df   = pd.read_csv(counts_path)
metadata_df = pd.read_csv(meta_path)

print("Counts loaded:", counts_df.shape)
print("Metadata loaded:", metadata_df.shape)

if "sample_id" in counts_df.columns:
    counts_df = counts_df.drop(columns=["sample_id"])

if "sample_id" not in metadata_df.columns:
    raise RuntimeError("metadata_df must contain sample_id column")

metadata_df.index = metadata_df["sample_id"].astype(str).values

adata = ad.AnnData(
    X   = counts_df.values.astype(np.float64),
    obs = metadata_df,
    var = pd.DataFrame(index=counts_df.columns)
)

print("AnnData created:")
print(adata)
print("X shape:", adata.X.shape)
print("obs columns:", list(adata.obs.columns))
')

cat("✓ AnnData created\n")

# =============================================================================
# Run scCODA model
# =============================================================================
cat("\nRunning scCODA model (HMC sampling)...\n")

py_run_string('
import numpy as np
import pandas as pd
from sccoda.util import comp_ana as ana

print("Running scCODA...")
print("adata X shape:", adata.X.shape)
print("adata obs columns:", list(adata.obs.columns))
print("adata var names (first 5):", list(adata.var_names)[:5])

if "condition_code" in adata.obs.columns:
    adata.obs["condition_code"] = pd.Categorical(
        adata.obs["condition_code"],
        categories=["HC", "UCSC", "UC"],
        ordered=True
    )

model = ana.CompositionalAnalysis(
    adata,
    formula="condition_code",
    reference_cell_type="automatic"
)

result = model.sample_hmc(
    num_results=20000,
    num_burnin=5000
)

print("scCODA finished.")

effect_df = result.effect_df
print("effect_df shape:", effect_df.shape)

print("Model summary:")
print(result.summary())
')

cat("✓ scCODA model complete\n")

# =============================================================================
# scCODA composition boxplots
# =============================================================================
cat("\nGenerating per-celltype composition plots...\n")

py_run_string('
import re, unicodedata, gc
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.close("all")
np.random.seed(42)

print("Counts path:", counts_path)
print("Meta path  :", meta_path)
print("Saving to  :", out_dir)

if "effect_df" not in globals():
    raise RuntimeError("effect_df not found. Run scCODA first.")

# Load counts + metadata
counts_raw = pd.read_csv(counts_path)
meta_raw   = pd.read_csv(meta_path)

if "sample_id" not in counts_raw.columns and "Unnamed: 0" in counts_raw.columns:
    counts_raw = counts_raw.rename(columns={"Unnamed: 0": "sample_id"})
if "sample_id" not in meta_raw.columns and "Unnamed: 0" in meta_raw.columns:
    meta_raw = meta_raw.rename(columns={"Unnamed: 0": "sample_id"})

if "sample_id" not in meta_raw.columns:
    raise RuntimeError("metadata must contain sample_id")
if "condition_code" not in meta_raw.columns:
    raise RuntimeError("metadata must contain condition_code")

meta_raw["sample_id"] = meta_raw["sample_id"].astype(str)

if "sample_id" in counts_raw.columns:
    counts_raw["sample_id"] = counts_raw["sample_id"].astype(str)
    df = counts_raw.merge(meta_raw[["sample_id", "condition_code"]], on="sample_id", how="left")
    if df["condition_code"].isna().any():
        missing = df.loc[df["condition_code"].isna(), "sample_id"].unique().tolist()
        raise RuntimeError(f"sample_id missing in metadata for: {missing}")
    cell_type_cols = [c for c in counts_raw.columns if c != "sample_id"]
else:
    if counts_raw.shape[0] != meta_raw.shape[0]:
        raise RuntimeError("Counts rows must match metadata rows if no sample_id column.")
    df = counts_raw.copy()
    df["sample_id"] = meta_raw["sample_id"].values
    df["condition_code"] = meta_raw["condition_code"].values
    cell_type_cols = counts_raw.columns.tolist()

cell_counts = df[cell_type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
totals = cell_counts.sum(axis=1).replace(0, np.nan)

proportions = (cell_counts.div(totals, axis=0) * 100.0).fillna(0.0)
proportions["condition"] = df["condition_code"].values
proportions["sample_id"] = df["sample_id"].values

cell_types = cell_type_cols

# Extract inclusion probabilities
def extract_ip(effect_df, covariate):
    x = effect_df.loc[covariate]
    if isinstance(x, pd.DataFrame) and "Inclusion probability" in x.columns:
        return x["Inclusion probability"].to_dict()
    if isinstance(x, pd.Series) and isinstance(x.index, pd.MultiIndex):
        for lvl in range(x.index.nlevels):
            vals = x.index.get_level_values(lvl).astype(str)
            if any(v == "Inclusion probability" for v in vals):
                ip = x.xs("Inclusion probability", level=lvl)
                return ip.to_dict()
    raise RuntimeError(f"Could not parse inclusion probabilities for {covariate}")

uc_ip_dict   = extract_ip(effect_df, "condition_code[T.UC]")
ucsc_ip_dict = extract_ip(effect_df, "condition_code[T.UCSC]")

# Match column names to scCODA keys
def norm_name(s):
    s = str(s).lower()
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s

uc_keys_norm   = {norm_name(k): k for k in uc_ip_dict.keys()}
ucsc_keys_norm = {norm_name(k): k for k in ucsc_ip_dict.keys()}

ct_to_uc_key, ct_to_ucsc_key = {}, {}
unmatched = []

for ct in cell_types:
    n = norm_name(ct)
    k1 = uc_keys_norm.get(n)
    k2 = ucsc_keys_norm.get(n)
    ct_to_uc_key[ct] = k1
    ct_to_ucsc_key[ct] = k2
    if k1 is None or k2 is None:
        unmatched.append(ct)

if unmatched:
    print("WARNING unmatched cell types (stars may be missing):")
    print(unmatched)

# Styling
FIG_W, FIG_H = 18, 13
DPI = 220

TITLE_FS      = 54
AXIS_LABEL_FS = 47
YTICK_FS      = 44
XTICK_FS      = 35
LEGEND_FS     = 22

Y_LABEL_PAD = 30
X_LABEL_PAD = 16

POINT_SIZE = 450
BOX_LW     = 3.5
WHISKER_LW = 3.2
CAP_LW     = 3.2
MEDIAN_LW  = 4.2

BRACKET_LW = 4.0
STAR_FS    = 44
GRID_ALPHA = 0.18

condition_order = ["HC", "UCSC", "UC"]
colors = {"HC": "#2166ac", "UCSC": "#92c5de", "UC": "#b2182b"}
GREEN = "#00BFA5"

def ip_to_stars(ip):
    if ip is None or (isinstance(ip, float) and (not np.isfinite(ip))):
        return None
    if ip >= 0.9:
        return "**"
    if ip >= 0.7:
        return "*"
    return None

def add_bracket(ax, x1, x2, y, h, color, text):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
            color=color, linewidth=BRACKET_LW, solid_capstyle="round")
    ax.text((x1+x2)/2, y+h + 0.12*h, text,
            ha="center", va="bottom", fontsize=STAR_FS, fontweight="bold", color=color)

_reserved = {
    "CON","PRN","AUX","NUL","COM1","COM2","COM3","COM4","COM5","COM6","COM7","COM8","COM9",
    "LPT1","LPT2","LPT3","LPT4","LPT5","LPT6","LPT7","LPT8","LPT9"
}

def safe_filename(name, maxlen=50):
    s = str(name)
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    s = s.replace("+", "plus")
    s = re.sub(r"\\s+", "_", s.strip())
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", s)
    s = s.strip(" ._-")
    if not s:
        s = "celltype"
    if s.upper() in _reserved:
        s = "X_" + s
    return s[:maxlen]

def q_mid(x, q):
    try:
        return float(np.quantile(x, q, method="midpoint"))
    except TypeError:
        return float(np.quantile(x, q, interpolation="midpoint"))

def midpoint_box_stats(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = np.sort(x)

    if x.size == 0:
        return dict(whislo=0.0, q1=0.0, med=0.0, q3=0.0, whishi=0.0, fliers=[])

    q1  = q_mid(x, 0.25)
    med = q_mid(x, 0.50)
    q3  = q_mid(x, 0.75)

    iqr = q3 - q1
    lo  = q1 - 1.5 * iqr
    hi  = q3 + 1.5 * iqr

    whislo = float(np.min(x[x >= lo])) if np.any(x >= lo) else float(np.min(x))
    whishi = float(np.max(x[x <= hi])) if np.any(x[x <= hi]) else float(np.max(x))

    return dict(whislo=whislo, q1=q1, med=med, q3=q3, whishi=whishi, fliers=[])

# Plot per cell type
made = 0
failed = []

for i, ct in enumerate(cell_types, start=1):
    fig, ax = plt.subplots(1, 1, figsize=(FIG_W, FIG_H))

    box_data   = []
    y_all      = []
    stats_list = []

    for cond in condition_order:
        vals = proportions.loc[proportions["condition"] == cond, ct].values.astype(float)
        box_data.append(vals)
        y_all.extend(vals.tolist())
        stats_list.append(midpoint_box_stats(vals))

    bxp = ax.bxp(
        stats_list,
        positions=[1, 2, 3],
        widths=0.72,
        showfliers=False,
        patch_artist=True,
        boxprops=dict(linewidth=BOX_LW, edgecolor="black"),
        whiskerprops=dict(linewidth=WHISKER_LW, color="black"),
        capprops=dict(linewidth=CAP_LW, color="black"),
        medianprops=dict(linewidth=MEDIAN_LW, color="black")
    )

    for patch, cond in zip(bxp["boxes"], condition_order):
        patch.set_facecolor(colors[cond])
        patch.set_alpha(0.90)

    for j, cond in enumerate(condition_order):
        vals = box_data[j]
        x_j = np.random.normal(loc=j+1, scale=0.06, size=len(vals))
        ax.scatter(x_j, vals, s=POINT_SIZE, c="black", alpha=0.95, zorder=5,
                   edgecolors="white", linewidths=1.6)

    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["HC", "UCSC", "UC"], fontsize=XTICK_FS, fontweight="bold")
    ax.tick_params(axis="y", labelsize=YTICK_FS)

    ax.set_ylabel("Proportion (%)", fontsize=AXIS_LABEL_FS, fontweight="bold", labelpad=Y_LABEL_PAD)
    ax.set_xlabel("Condition", fontsize=AXIS_LABEL_FS, fontweight="bold", labelpad=X_LABEL_PAD)

    ax.yaxis.grid(True, linestyle="--", alpha=GRID_ALPHA)
    ax.set_axisbelow(True)

    ax.set_title(str(ct), fontsize=TITLE_FS, fontweight="bold", pad=26)

    y_max = max(y_all) if y_all else 1.0
    y_min = min(y_all) if y_all else 0.0
    span  = max(1.0, y_max - y_min)
    ax.set_ylim(0, y_max + 1.8 * span)

    for spine in ax.spines.values():
        spine.set_linewidth(2.8)

    uc_key   = ct_to_uc_key.get(ct)
    ucsc_key = ct_to_ucsc_key.get(ct)

    uc_ip   = float(uc_ip_dict.get(uc_key, np.nan)) if uc_key else np.nan
    ucsc_ip = float(ucsc_ip_dict.get(ucsc_key, np.nan)) if ucsc_key else np.nan

    base_y = y_max + 0.10 * span
    h      = 0.12 * span

    s_uc = ip_to_stars(uc_ip)
    if s_uc:
        add_bracket(ax, 1, 3, base_y, h, color=colors["UC"], text=s_uc)

    s_ucsc = ip_to_stars(ucsc_ip)
    if s_ucsc:
        add_bracket(ax, 1, 2, base_y + 2.8 * h, h, color=GREEN, text=s_ucsc)

    legend_elements = [
        Patch(facecolor=colors["HC"],   edgecolor="black", alpha=0.90, label="Healthy Control (HC)"),
        Patch(facecolor=colors["UCSC"], edgecolor="black", alpha=0.90, label="UC Self-Control (UCSC)"),
        Patch(facecolor=colors["UC"],   edgecolor="black", alpha=0.90, label="Ulcerative Colitis (UC)"),
        Line2D([0], [0], color=colors["UC"], linewidth=6, label="scCODA: UC vs HC"),
        Line2D([0], [0], color=GREEN,        linewidth=6, label="scCODA: UCSC vs HC"),
        Line2D([0], [0], color="white", marker="None", label="Stars show IP:  * >= 0.7   ** >= 0.9"),
    ]

    leg = ax.legend(
        handles=legend_elements,
        loc="upper right",
        fontsize=LEGEND_FS,
        frameon=True,
        fancybox=True,
        framealpha=0.96,
        borderpad=0.5,
        labelspacing=0.35,
        handletextpad=0.5,
        handlelength=2.6,
        handleheight=0.9
    )
    leg.get_frame().set_linewidth(2.0)

    fig.tight_layout()

    fname = f"{i:02d}_sccoda_{safe_filename(ct)}.png"
    out_path = out_dir / fname

    try:
        fig.savefig(str(out_path), dpi=DPI, bbox_inches="tight", facecolor="white")
        made += 1
    except Exception as e:
        failed.append((ct, str(e)))
        print(f"ERROR saving {ct}: {e}")

    plt.close(fig)
    gc.collect()

print(f"Done. Saved {made} plots to: {out_dir}")

if failed:
    print("FAILED plots:")
    for ct, err in failed:
        print(f"  {ct} -> {err}")
')

cat("✓ scCODA plots complete\n")
