#!/usr/bin/env python3
"""
multiomics_ml.py
----------------
Integrates four genomic data layers into a unified bin-level feature matrix,
then trains machine learning models to correlate methylation + 3D genome
organisation (Hi-C) with gene expression.

Input files
-----------
  1. --annot-methyl   BED/TSV from annotate_methylation_gff.py
                      Per-feature methylation (gene-level resolution).
                      Columns: replicon, feature_type, feat_start, feat_end,
                               strand, feat_length_bp, locus_tag, gene, product,
                               mod_code, mod_label, mod_description,
                               raw_count, norm_count

  2. --bin-methyl     BED from methylation_windows.py
                      Methylation counts in fixed genomic bins.
                      Columns: replicon, bin_start, bin_end, strand,
                               mod_code, mod_label, raw_count, norm_count

  3. --deseq2         TSV from DESeq2
                      Differential expression per gene.
                      Columns: gene_id, baseMean, log2FoldChange, lfcSE,
                               stat, pvalue, padj

  4. --hic            BED from a .cool file
                      Hi-C contact values per bin.
                      Columns: chrom, start, end, contact_value

Harmonisation strategy
----------------------
  A master bin grid is built from the Hi-C bins (or at --resolution if
  provided).  For each master bin:

    Methylation features  (from --bin-methyl):
      For every (mod_label × strand) combination:
        • methyl_<mod>_<strand>_raw_sum   — total raw methylation count
        • methyl_<mod>_<strand>_norm_mean — mean normalised methylation
      Aggregated from all input bins whose midpoint falls in the master bin.

    Hi-C features  (from --hic):
      • hic_contact_mean   — mean contact value across overlapping Hi-C bins
      • hic_contact_max    — max  contact value
      • hic_contact_sum    — sum  contact value

    Expression features  (from --deseq2 + --annot-methyl for gene positions):
      Genes are assigned to master bins by midpoint of [feat_start, feat_end).
        • expr_log2fc_mean    — mean log2FoldChange of genes in the bin
        • expr_log2fc_max     — max  |log2FoldChange|
        • expr_baseMean_mean  — mean baseMean
        • expr_padj_min       — min padj (most significant gene in bin)
        • expr_n_genes        — number of DEGs in the bin
        • expr_n_sig_genes    — number with padj < 0.05

    Gene-level annotation methylation features (from --annot-methyl):
      For every (mod_label) combination:
        • annot_<mod>_raw_sum
        • annot_<mod>_norm_mean

    Spatial / positional features:
      • replicon_enc  — integer-encoded replicon
      • bin_center    — genomic centre of the bin (bp)

Target variable
---------------
  expr_log2fc_mean  (mean log2FoldChange of genes in the bin)
  Only bins that contain at least one gene with expression data are used.

Models trained
--------------
  1. XGBoost (gradient boosting) — primary model with SHAP
  2. Random Forest               — comparison model
  Both are evaluated on a held-out test set (--test-frac, default 0.2).

Outputs (saved to --output-dir)
-------------------------------
  master_bins.tsv           — full harmonised feature matrix
  model_metrics.tsv         — RMSE, R², Pearson r for both models
  shap_summary.tsv          — mean |SHAP| per feature
  shap_values.tsv           — full SHAP matrix (bins × features)
  rf_importance.tsv         — Random Forest feature importances
  plots/
    shap_importance.pdf/png     — top-N SHAP bar chart
    shap_beeswarm.pdf/png       — SHAP summary (beeswarm)
    shap_dependence_*.pdf/png   — dependence plots for top features
    shap_heatmap.pdf/png        — bin × feature SHAP heatmap
    rf_importance.pdf/png       — RF importance bar chart
    hic_vs_expr.pdf/png         — scatter: Hi-C contact vs expression
    methyl_vs_expr.pdf/png      — scatter: methylation vs expression
    model_comparison.pdf/png    — XGB vs RF test predictions

Usage
-----
  python multiomics_ml.py \\
      --annot-methyl annot_methyl.tsv \\
      --bin-methyl   bin_methyl_6mA.bed \\
      --deseq2       deseq2_results.tsv \\
      --hic          hic_contacts.bed \\
      --output-dir   ml_results/

  # Multi-file methylation (one per mod type, will be concatenated)
  python multiomics_ml.py \\
      --annot-methyl annot_methyl_6mA.tsv annot_methyl_5mC.tsv \\
      --bin-methyl   bin_methyl_6mA.bed bin_methyl_5mC.bed \\
      --deseq2       deseq2_results.tsv \\
      --hic          hic_contacts.bed \\
      --output-dir   ml_results/

  # Override master bin resolution (bp)
  python multiomics_ml.py ... --resolution 10000

Dependencies
------------
  numpy, pandas, scipy, scikit-learn, xgboost, shap, matplotlib, seaborn
  Install:
    pip install numpy pandas scipy scikit-learn xgboost shap matplotlib seaborn
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ── Optional imports ──────────────────────────────────────────────────────────
try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import LabelEncoder

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ===========================================================================
# 1.  Loaders
# ===========================================================================

def _read_flexible(paths: list[Path], comment: str = "#") -> pd.DataFrame:
    """
    Read one or more TSV/BED files and concatenate them.

    Handles two header styles produced by the pipeline tools:
      - Normal header:   replicon\tbin_start\t...
      - Commented header: # replicon\tbin_start\t...   (methylation_windows.py style)

    For commented headers the '#' is stripped and the line is used as the
    column header; subsequent lines starting with '#' are treated as comments.
    """
    frames = []
    for p in paths:
        # Peek at the first non-empty line to detect a commented header
        header_line = None
        with open(p, "r", encoding="utf-8") as fh:
            for raw in fh:
                stripped = raw.strip()
                if stripped:
                    header_line = stripped
                    break

        if header_line is None:
            continue  # empty file

        if header_line.startswith("#"):
            # The header itself is commented — strip the '#', read manually
            col_names = header_line.lstrip("#").strip().split("\t")
            df = pd.read_csv(
                p, sep="\t", comment="#",
                names=col_names, header=None,
                dtype=str, low_memory=False,
                skiprows=lambda i: False,   # we handle skipping via comment
            )
            # The first data row might be the header line itself if pandas
            # didn't skip it (depends on whether pandas treated it as comment).
            # Drop any row whose first cell equals the first column name.
            if len(df) > 0 and df.iloc[0, 0] == col_names[0]:
                df = df.iloc[1:].reset_index(drop=True)
        else:
            # Normal header — read directly, skip lines starting with '#'
            df = pd.read_csv(
                p, sep="\t", comment="#",
                dtype=str, low_memory=False,
            )
            df.columns = [c.lstrip("#").strip() for c in df.columns]

        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def load_annot_methyl(paths: list[Path]) -> pd.DataFrame:
    """Load gene-level annotation methylation files."""
    df = _read_flexible(paths)
    required = {"replicon", "feat_start", "feat_end", "strand",
                "mod_label", "raw_count", "norm_count", "locus_tag"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"annot-methyl missing columns: {missing}")
    df["feat_start"] = pd.to_numeric(df["feat_start"], errors="coerce")
    df["feat_end"]   = pd.to_numeric(df["feat_end"],   errors="coerce")
    df["raw_count"]  = pd.to_numeric(df["raw_count"],  errors="coerce").fillna(0)
    df["norm_count"] = pd.to_numeric(df["norm_count"], errors="coerce").fillna(0)
    df["feat_mid"]   = (df["feat_start"] + df["feat_end"]) / 2
    logger.info("annot-methyl: %d rows, %d genes, mod types: %s",
                len(df), df["locus_tag"].nunique(),
                sorted(df["mod_label"].dropna().unique()))
    return df


def load_bin_methyl(paths: list[Path]) -> pd.DataFrame:
    """Load binned methylation files."""
    df = _read_flexible(paths)
    required = {"replicon", "bin_start", "bin_end", "strand",
                "mod_label", "raw_count", "norm_count"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"bin-methyl missing columns: {missing}")
    df["bin_start"] = pd.to_numeric(df["bin_start"], errors="coerce")
    df["bin_end"]   = pd.to_numeric(df["bin_end"],   errors="coerce")
    df["raw_count"] = pd.to_numeric(df["raw_count"], errors="coerce").fillna(0)
    df["norm_count"]= pd.to_numeric(df["norm_count"],errors="coerce").fillna(0)
    df["bin_mid"]   = (df["bin_start"] + df["bin_end"]) / 2
    logger.info("bin-methyl: %d rows, mod types: %s, strands: %s",
                len(df),
                sorted(df["mod_label"].dropna().unique()),
                sorted(df["strand"].dropna().unique()))
    return df


def load_deseq2(path: Path) -> pd.DataFrame:
    """Load DESeq2 results. Flexible column name matching."""
    df = _read_flexible([path])
    # Rename common variants
    rename_map = {
        "gene": "gene_id", "GeneID": "gene_id", "ID": "gene_id",
        "locus_tag": "gene_id",
        "log2FC": "log2FoldChange", "LFC": "log2FoldChange",
        "basemean": "baseMean",
        "adj.P.Val": "padj", "FDR": "padj",
    }
    df.rename(columns={k: v for k, v in rename_map.items()
                        if k in df.columns and v not in df.columns},
              inplace=True)
    required = {"gene_id", "log2FoldChange", "baseMean"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DESeq2 missing columns: {missing}")
    for col in ["log2FoldChange", "baseMean", "padj", "pvalue"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["log2FoldChange"])
    logger.info("DESeq2: %d genes, LFC range [%.3f, %.3f]",
                len(df), df["log2FoldChange"].min(), df["log2FoldChange"].max())
    return df


def load_hic(path: Path) -> pd.DataFrame:
    """Load Hi-C contact BED file."""
    df = _read_flexible([path])
    rename_map = {"chr": "chrom", "chromosome": "chrom",
                  "value": "contact_value", "score": "contact_value"}
    df.rename(columns={k: v for k, v in rename_map.items()
                        if k in df.columns and v not in df.columns},
              inplace=True)
    required = {"chrom", "start", "end", "contact_value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"hic missing columns: {missing}")
    df["start"]         = pd.to_numeric(df["start"],         errors="coerce")
    df["end"]           = pd.to_numeric(df["end"],           errors="coerce")
    df["contact_value"] = pd.to_numeric(df["contact_value"], errors="coerce").fillna(0)
    df["hic_mid"]       = (df["start"] + df["end"]) / 2
    logger.info("Hi-C: %d bins, %d replicon(s), contact range [%.2f, %.2f]",
                len(df), df["chrom"].nunique(),
                df["contact_value"].min(), df["contact_value"].max())
    return df


# ===========================================================================
# 2.  Master bin grid
# ===========================================================================

def build_master_bins(
    hic_df: pd.DataFrame,
    resolution: int | None,
) -> pd.DataFrame:
    """
    Build the master bin grid.

    If --resolution is given, bins of that size are created from 0 to
    max(end) for each replicon present in the Hi-C data.
    Otherwise the Hi-C bins themselves are used as the master grid.

    Returns a DataFrame with columns: replicon, bin_start, bin_end, bin_center.
    """
    if resolution is not None:
        rows = []
        for chrom, grp in hic_df.groupby("chrom"):
            max_end = int(grp["end"].max())
            for s in range(0, max_end, resolution):
                rows.append({
                    "replicon":   chrom,
                    "bin_start":  s,
                    "bin_end":    min(s + resolution, max_end),
                })
        master = pd.DataFrame(rows)
    else:
        master = hic_df[["chrom", "start", "end"]].copy()
        master.columns = ["replicon", "bin_start", "bin_end"]
        master = master.drop_duplicates().reset_index(drop=True)

    master["bin_center"] = (master["bin_start"] + master["bin_end"]) / 2
    master = master.sort_values(["replicon", "bin_start"]).reset_index(drop=True)
    logger.info("Master grid: %d bins across %d replicon(s)",
                len(master), master["replicon"].nunique())
    return master


# ===========================================================================
# 3.  Aggregation helpers
# ===========================================================================

def _assign_to_master(
    data_df: pd.DataFrame,
    master: pd.DataFrame,
    data_replicon_col: str,
    data_pos_col: str,
) -> pd.Series:
    """
    For each row in data_df, find the master bin index whose
    [bin_start, bin_end) contains data_pos_col on the same replicon.

    Returns a Series of master bin indices (NaN if no matching bin).
    Uses a vectorised merge-asof approach per replicon.
    """
    result = pd.Series(np.nan, index=data_df.index, dtype=float)

    for rep, grp in data_df.groupby(data_replicon_col):
        m = master[master["replicon"] == rep].reset_index()
        if m.empty:
            continue
        # sorted master bins
        m_sorted = m.sort_values("bin_start")

        positions = grp[data_pos_col].values
        bin_starts = m_sorted["bin_start"].values
        bin_ends   = m_sorted["bin_end"].values
        orig_idx   = m_sorted["index"].values   # master DataFrame index

        # Binary search: leftmost bin_start <= position
        insert = np.searchsorted(bin_starts, positions, side="right") - 1
        valid  = (insert >= 0)

        for i, (pos, bi, v) in enumerate(zip(positions, insert, valid)):
            if not v:
                continue
            if pos < bin_ends[bi]:   # position is inside the bin
                result.iloc[grp.index[i]] = orig_idx[bi]

    return result


def aggregate_bin_methyl(
    bin_methyl: pd.DataFrame,
    master: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate methylation bins -> master bins.

    For each (mod_label × strand) combination:
        methyl_<mod>_<strand>_raw_sum
        methyl_<mod>_<strand>_norm_mean
    """
    bin_methyl = bin_methyl.copy()
    bin_methyl["_master_idx"] = _assign_to_master(
        bin_methyl, master, "replicon", "bin_mid"
    )
    bin_methyl = bin_methyl.dropna(subset=["_master_idx"])
    bin_methyl["_master_idx"] = bin_methyl["_master_idx"].astype(int)

    # Sanitise label for column names
    bin_methyl["_col_key"] = (
        "methyl_"
        + bin_methyl["mod_label"].str.replace(r"[^A-Za-z0-9]", "", regex=True)
        + "_"
        + bin_methyl["strand"].str.replace(r"[^A-Za-z0-9]", "dot", regex=True)
    )

    agg_rows = []
    for master_idx, grp in bin_methyl.groupby("_master_idx"):
        row = {"_master_idx": master_idx}
        for key, sub in grp.groupby("_col_key"):
            row[f"{key}_raw_sum"]   = sub["raw_count"].sum()
            row[f"{key}_norm_mean"] = sub["norm_count"].mean()
        agg_rows.append(row)

    if not agg_rows:
        return master[["replicon", "bin_start", "bin_end", "bin_center"]].copy()

    agg = pd.DataFrame(agg_rows).fillna(0)
    result = master.rename_axis("_master_idx").reset_index().merge(agg, on="_master_idx", how="left")
    methyl_cols = [c for c in result.columns if c.startswith("methyl_")]
    result[methyl_cols] = result[methyl_cols].fillna(0)
    logger.info("Methylation aggregation: %d feature column(s)", len(methyl_cols))
    return result.drop(columns=["_master_idx", "index"], errors="ignore")


def aggregate_hic(
    hic_df: pd.DataFrame,
    master: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate Hi-C contact values -> master bins.
    Produces: hic_contact_mean, hic_contact_max, hic_contact_sum.
    """
    hic_df = hic_df.copy()
    hic_df["_master_idx"] = _assign_to_master(
        hic_df, master, "chrom", "hic_mid"
    )
    hic_df = hic_df.dropna(subset=["_master_idx"])
    hic_df["_master_idx"] = hic_df["_master_idx"].astype(int)

    agg = (hic_df.groupby("_master_idx")["contact_value"]
           .agg(hic_contact_mean="mean",
                hic_contact_max="max",
                hic_contact_sum="sum")
           .reset_index())

    result = master.rename_axis("_master_idx").reset_index().merge(agg, on="_master_idx", how="left")
    for col in ["hic_contact_mean", "hic_contact_max", "hic_contact_sum"]:
        result[col] = result[col].fillna(0)
    logger.info("Hi-C aggregation: %d bins with contact data",
                agg["_master_idx"].nunique())
    return result.drop(columns=["_master_idx", "index"], errors="ignore")


def aggregate_expression(
    annot_methyl: pd.DataFrame,
    deseq2: pd.DataFrame,
    master: pd.DataFrame,
) -> pd.DataFrame:
    """
    Map genes to master bins using gene midpoint, then aggregate expression.

    Genes are located via the annot-methyl file (which contains coordinates).
    Expression values come from the DESeq2 file, joined on locus_tag / gene_id.

    Produces per-bin:
        expr_log2fc_mean, expr_log2fc_max_abs, expr_baseMean_mean,
        expr_padj_min, expr_n_genes, expr_n_sig_genes
    """
    # Keep unique gene positions from annot_methyl
    gene_pos = (annot_methyl[["locus_tag", "replicon", "feat_mid"]]
                .drop_duplicates(subset=["locus_tag"]))

    # Join expression
    expr = deseq2.rename(columns={"gene_id": "locus_tag"})
    gene_expr = gene_pos.merge(expr, on="locus_tag", how="inner")

    if gene_expr.empty:
        logger.warning("No gene IDs matched between annot-methyl and DESeq2.")
        return master[["replicon", "bin_start", "bin_end", "bin_center"]].copy()

    gene_expr["_master_idx"] = _assign_to_master(
        gene_expr, master, "replicon", "feat_mid"
    )
    gene_expr = gene_expr.dropna(subset=["_master_idx"])
    gene_expr["_master_idx"] = gene_expr["_master_idx"].astype(int)

    def sig_count(x):
        return (x < 0.05).sum()

    agg = (gene_expr.groupby("_master_idx")
           .agg(
               expr_log2fc_mean   =("log2FoldChange", "mean"),
               expr_log2fc_max_abs=("log2FoldChange", lambda x: x.abs().max()),
               expr_baseMean_mean =("baseMean",        "mean"),
               expr_padj_min      =("padj",            "min"),
               expr_n_genes       =("locus_tag",       "count"),
               expr_n_sig_genes   =("padj",            sig_count),
           )
           .reset_index())

    result = master.rename_axis("_master_idx").reset_index().merge(agg, on="_master_idx", how="left")
    logger.info("Expression aggregation: %d bins contain >= 1 gene",
                agg["_master_idx"].nunique())
    return result.drop(columns=["_master_idx", "index"], errors="ignore")


def aggregate_annot_methyl(
    annot_methyl: pd.DataFrame,
    master: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate gene-level annotation methylation -> master bins.
    Produces: annot_<mod>_raw_sum, annot_<mod>_norm_mean per mod_label.
    """
    am = annot_methyl.copy()
    am["_master_idx"] = _assign_to_master(am, master, "replicon", "feat_mid")
    am = am.dropna(subset=["_master_idx"])
    am["_master_idx"] = am["_master_idx"].astype(int)

    am["_col_key"] = "annot_" + am["mod_label"].str.replace(
        r"[^A-Za-z0-9]", "", regex=True
    )

    agg_rows = []
    for master_idx, grp in am.groupby("_master_idx"):
        row = {"_master_idx": master_idx}
        for key, sub in grp.groupby("_col_key"):
            row[f"{key}_raw_sum"]   = sub["raw_count"].sum()
            row[f"{key}_norm_mean"] = sub["norm_count"].mean()
        agg_rows.append(row)

    if not agg_rows:
        return master[["replicon", "bin_start", "bin_end", "bin_center"]].copy()

    agg = pd.DataFrame(agg_rows).fillna(0)
    result = master.rename_axis("_master_idx").reset_index().merge(agg, on="_master_idx", how="left")
    annot_cols = [c for c in result.columns if c.startswith("annot_")]
    result[annot_cols] = result[annot_cols].fillna(0)
    logger.info("Annot-methyl aggregation: %d feature column(s)", len(annot_cols))
    return result.drop(columns=["_master_idx", "index"], errors="ignore")


# ===========================================================================
# 4.  Build master feature matrix
# ===========================================================================

def build_feature_matrix(
    master: pd.DataFrame,
    bin_methyl_agg: pd.DataFrame,
    hic_agg: pd.DataFrame,
    expr_agg: pd.DataFrame,
    annot_methyl_agg: pd.DataFrame,
) -> pd.DataFrame:
    """
    Left-join all aggregated layers onto the master bin grid.
    Adds replicon_enc (integer-encoded replicon).
    Returns the full feature matrix (including target column expr_log2fc_mean).
    """
    key = ["replicon", "bin_start", "bin_end", "bin_center"]

    df = master[key].copy()

    def _merge(left, right, suffix):
        right_cols = [c for c in right.columns if c not in key or c in key[:3]]
        return left.merge(
            right[[c for c in right.columns if c in key or c not in left.columns]],
            on=["replicon", "bin_start", "bin_end"],
            how="left",
        )

    df = df.merge(
        bin_methyl_agg.drop(columns=["bin_center"], errors="ignore"),
        on=["replicon", "bin_start", "bin_end"], how="left",
    )
    df = df.merge(
        hic_agg.drop(columns=["bin_center"], errors="ignore"),
        on=["replicon", "bin_start", "bin_end"], how="left",
    )
    df = df.merge(
        expr_agg.drop(columns=["bin_center"], errors="ignore"),
        on=["replicon", "bin_start", "bin_end"], how="left",
    )
    df = df.merge(
        annot_methyl_agg.drop(columns=["bin_center"], errors="ignore"),
        on=["replicon", "bin_start", "bin_end"], how="left",
    )

    # Encode replicon
    le = LabelEncoder()
    df["replicon_enc"] = le.fit_transform(df["replicon"].fillna("unknown"))

    # Fill numeric NaN with 0 (absent methylation / contact = 0)
    num_cols = df.select_dtypes(include=[np.number]).columns
    df[num_cols] = df[num_cols].fillna(0)

    logger.info("Feature matrix: %d bins × %d columns", len(df), len(df.columns))
    return df


# ===========================================================================
# 5.  ML models
# ===========================================================================

TARGET = "expr_log2fc_mean"

EXCLUDE_COLS = {
    "replicon", "bin_start", "bin_end", "bin_center",
    # expression columns that would be data leakage
    "expr_log2fc_mean", "expr_log2fc_max_abs",
    "expr_padj_min", "expr_baseMean_mean",
    "expr_n_genes", "expr_n_sig_genes",
}


def prepare_ml_data(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series,
                                                 list[str]]:
    """
    Subset to bins with expression data, select feature columns, drop NaN target.
    Returns X (features), y (target), feature_names.
    """
    model_df = df.dropna(subset=[TARGET]).copy()
    model_df = model_df[model_df["expr_n_genes"] > 0]

    if len(model_df) < 4:
        raise ValueError(
            f"Only {len(model_df)} bin(s) have expression data (need >= 4). "
            "Check that gene IDs in DESeq2 match locus_tag values in annot-methyl."
        )

    feature_cols = [
        c for c in model_df.columns
        if c not in EXCLUDE_COLS
        and model_df[c].dtype in [np.float64, np.int64, np.float32, np.int32]
    ]

    X = model_df[feature_cols].fillna(0)
    y = model_df[TARGET]

    logger.info("ML dataset: %d samples × %d features", len(X), len(X.columns))
    return X, y, feature_cols


def train_xgboost(
    X_train, y_train, X_test, y_test, seed: int
) -> tuple[object, dict]:
    if not HAS_XGB:
        raise ImportError("xgboost is not installed.")

    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest  = xgb.DMatrix(X_test,  label=y_test)

    params = dict(
        booster="gbtree", objective="reg:squarederror",
        eta=0.05, max_depth=4, subsample=0.8,
        colsample_bytree=0.8, min_child_weight=3,
        gamma=0.1, lambda_=1.0, alpha=0.1, seed=seed,
    )
    model = xgb.train(
        params, dtrain, num_boost_round=500,
        evals=[(dtest, "test")],
        early_stopping_rounds=30,
        verbose_eval=False,
    )
    pred = model.predict(xgb.DMatrix(X_test))
    metrics = _metrics(y_test.values, pred, "XGBoost")
    return model, metrics


def train_rf(
    X_train, y_train, X_test, y_test, seed: int
) -> tuple[RandomForestRegressor, dict]:
    model = RandomForestRegressor(
        n_estimators=500, max_depth=None,
        min_samples_leaf=3, n_jobs=-1, random_state=seed,
    )
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    metrics = _metrics(y_test.values, pred, "RandomForest")
    return model, metrics


def _metrics(y_true, y_pred, name: str) -> dict:
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    r2   = float(r2_score(y_true, y_pred))
    r, _ = pearsonr(y_true, y_pred)
    logger.info("%s — RMSE: %.4f  R²: %.4f  Pearson r: %.4f", name, rmse, r2, r)
    return {"model": name, "rmse": rmse, "r_squared": r2, "pearson_r": float(r)}


def compute_shap(model, X_train: pd.DataFrame) -> np.ndarray:
    """Compute SHAP values for an XGBoost model."""
    if not HAS_SHAP:
        logger.warning("shap not installed — skipping SHAP computation.")
        return None
    explainer = shap.TreeExplainer(model)
    shap_vals = explainer.shap_values(X_train)
    return shap_vals


# ===========================================================================
# 6.  Save outputs
# ===========================================================================

def save_outputs(
    out_dir: Path,
    master_df: pd.DataFrame,
    X: pd.DataFrame,
    y: pd.Series,
    xgb_model,
    rf_model: RandomForestRegressor,
    xgb_metrics: dict,
    rf_metrics: dict,
    shap_vals: np.ndarray | None,
    fmt: str,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    # ── Master feature matrix ────────────────────────────────────────────────
    master_df.to_csv(out_dir / "master_bins.tsv", sep="\t", index=False)
    logger.info("Saved: master_bins.tsv")

    # ── Model metrics ────────────────────────────────────────────────────────
    metrics_df = pd.DataFrame([xgb_metrics, rf_metrics])
    metrics_df.to_csv(out_dir / "model_metrics.tsv", sep="\t", index=False)
    logger.info("Saved: model_metrics.tsv")

    # ── RF importance ────────────────────────────────────────────────────────
    rf_imp = pd.DataFrame({
        "feature":    X.columns,
        "importance": rf_model.feature_importances_,
    }).sort_values("importance", ascending=False)
    rf_imp.to_csv(out_dir / "rf_importance.tsv", sep="\t", index=False)
    logger.info("Saved: rf_importance.tsv")

    # ── SHAP ─────────────────────────────────────────────────────────────────
    if shap_vals is not None:
        mean_abs = np.abs(shap_vals).mean(axis=0)
        shap_summary = pd.DataFrame({
            "feature":       X.columns,
            "mean_abs_shap": mean_abs,
        }).sort_values("mean_abs_shap", ascending=False)
        shap_summary.to_csv(out_dir / "shap_summary.tsv", sep="\t", index=False)

        shap_df = pd.DataFrame(shap_vals, columns=X.columns)
        shap_df.to_csv(out_dir / "shap_values.tsv", sep="\t", index=False)
        logger.info("Saved: shap_summary.tsv, shap_values.tsv")
    else:
        shap_summary = None

    if not HAS_PLOT:
        logger.warning("matplotlib/seaborn not installed — skipping plots.")
        return

    _make_plots(plot_dir, X, y, xgb_model, rf_model,
                shap_vals, shap_summary, rf_imp, master_df, fmt)


def _make_plots(
    plot_dir: Path,
    X: pd.DataFrame,
    y: pd.Series,
    xgb_model,
    rf_model: RandomForestRegressor,
    shap_vals: np.ndarray | None,
    shap_summary: pd.DataFrame | None,
    rf_imp: pd.DataFrame,
    master_df: pd.DataFrame,
    fmt: str,
    top_n: int = 20,
) -> None:

    def save(name: str):
        path = plot_dir / f"{name}.{fmt}"
        plt.tight_layout()
        plt.savefig(path, dpi=150, bbox_inches="tight")
        plt.close()
        logger.info("  Plot: %s", path.name)

    sns.set_style("whitegrid")
    pal = sns.color_palette("tab10")

    # ── 1. RF importance ─────────────────────────────────────────────────────
    top_rf = rf_imp.head(top_n)
    fig, ax = plt.subplots(figsize=(10, max(5, top_n * 0.4)))
    ax.barh(top_rf["feature"][::-1], top_rf["importance"][::-1], color=pal[1])
    ax.set_xlabel("Mean decrease in impurity")
    ax.set_title(f"Random Forest Feature Importance (top {top_n})")
    save("rf_importance")

    # ── 2. SHAP importance bar ────────────────────────────────────────────────
    if shap_summary is not None:
        top_shap = shap_summary.head(top_n)
        fig, ax = plt.subplots(figsize=(10, max(5, top_n * 0.4)))
        ax.barh(top_shap["feature"][::-1], top_shap["mean_abs_shap"][::-1],
                color=pal[0])
        ax.set_xlabel("Mean |SHAP value|")
        ax.set_title(f"XGBoost SHAP Importance (top {top_n})")
        save("shap_importance")

    # ── 3. SHAP beeswarm ─────────────────────────────────────────────────────
    if shap_vals is not None and HAS_SHAP:
        top_feat_idx = np.argsort(np.abs(shap_vals).mean(axis=0))[::-1][:top_n]
        X_top  = X.iloc[:, top_feat_idx]
        sv_top = shap_vals[:, top_feat_idx]
        fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.45)))
        shap.summary_plot(sv_top, X_top, plot_type="dot",
                          show=False, color_bar=True)
        plt.title("SHAP Summary (top features)")
        save("shap_beeswarm")

    # ── 4. SHAP heatmap ───────────────────────────────────────────────────────
    if shap_vals is not None:
        top_n_heat = min(10, shap_vals.shape[1])
        top_idx  = np.argsort(np.abs(shap_vals).mean(axis=0))[::-1][:top_n_heat]
        sv_heat  = shap_vals[:, top_idx]
        feat_names = X.columns[top_idx]
        heat_df  = pd.DataFrame(sv_heat, columns=feat_names)
        fig, ax  = plt.subplots(figsize=(max(8, len(heat_df) * 0.2), 5))
        sns.heatmap(heat_df.T, center=0, cmap="RdBu_r",
                    xticklabels=False, ax=ax,
                    cbar_kws={"label": "SHAP value"})
        ax.set_title("SHAP Heatmap (top features × bins)")
        ax.set_xlabel("Bins")
        save("shap_heatmap")

    # ── 5. SHAP dependence for top 3 features ────────────────────────────────
    if shap_vals is not None and HAS_SHAP:
        top3 = np.argsort(np.abs(shap_vals).mean(axis=0))[::-1][:3]
        for fi in top3:
            fname = X.columns[fi]
            fig, ax = plt.subplots(figsize=(7, 5))
            shap.dependence_plot(fi, shap_vals, X, show=False, ax=ax)
            ax.set_title(f"SHAP Dependence: {fname}")
            safe = fname.replace("/", "_").replace(" ", "_")
            save(f"shap_dependence_{safe}")

    # ── 6. Hi-C contact vs expression ────────────────────────────────────────
    plot_df = master_df.dropna(subset=["expr_log2fc_mean", "hic_contact_mean"])
    if not plot_df.empty:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(plot_df["hic_contact_mean"], plot_df["expr_log2fc_mean"],
                   alpha=0.5, s=15, color=pal[2])
        # regression line
        m, b = np.polyfit(plot_df["hic_contact_mean"],
                          plot_df["expr_log2fc_mean"], 1)
        xr = np.linspace(plot_df["hic_contact_mean"].min(),
                          plot_df["hic_contact_mean"].max(), 100)
        ax.plot(xr, m * xr + b, color="black", linewidth=1, linestyle="--")
        ax.set_xlabel("Mean Hi-C contact value (per bin)")
        ax.set_ylabel("Mean log2 Fold Change (per bin)")
        ax.set_title("Hi-C Contact vs Gene Expression")
        save("hic_vs_expr")

    # ── 7. Methylation vs expression ─────────────────────────────────────────
    methyl_norm_cols = [c for c in master_df.columns if "norm_mean" in c
                        and c.startswith("methyl_")]
    if methyl_norm_cols:
        # Use the first mod type found
        mc = methyl_norm_cols[0]
        plot_df2 = master_df.dropna(subset=["expr_log2fc_mean", mc])
        if not plot_df2.empty:
            fig, ax = plt.subplots(figsize=(7, 5))
            sc = ax.scatter(plot_df2[mc], plot_df2["expr_log2fc_mean"],
                            c=plot_df2.get("hic_contact_mean",
                                           pd.Series(0, index=plot_df2.index)),
                            cmap="viridis", alpha=0.6, s=15)
            plt.colorbar(sc, ax=ax, label="Hi-C contact (mean)")
            ax.set_xlabel(f"Norm methylation mean ({mc})")
            ax.set_ylabel("Mean log2 Fold Change")
            ax.set_title("Methylation vs Expression (coloured by Hi-C)")
            save("methyl_vs_expr")

    # ── 8. Model comparison: XGBoost vs RF predictions ───────────────────────
    if HAS_XGB and xgb_model is not None:
        xgb_pred = xgb_model.predict(xgb.DMatrix(X))
        rf_pred  = rf_model.predict(X)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, pred, name, colour in zip(
            axes,
            [xgb_pred, rf_pred],
            ["XGBoost", "Random Forest"],
            [pal[0], pal[1]],
        ):
            ax.scatter(y, pred, alpha=0.5, s=12, color=colour)
            lim = [min(y.min(), pred.min()), max(y.max(), pred.max())]
            ax.plot(lim, lim, "k--", linewidth=0.8)
            r2v = r2_score(y, pred)
            ax.set_title(f"{name}  (R²={r2v:.3f})")
            ax.set_xlabel("Observed log2FC")
            ax.set_ylabel("Predicted log2FC")
        fig.suptitle("Model Predictions vs Observed", fontsize=13)
        save("model_comparison")


# ===========================================================================
# 7.  Summary
# ===========================================================================

def print_summary(master_df: pd.DataFrame, xgb_m: dict, rf_m: dict) -> None:
    sep  = "=" * 66
    sep2 = "-" * 66
    methyl_cols = [c for c in master_df.columns if c.startswith("methyl_")]
    annot_cols  = [c for c in master_df.columns if c.startswith("annot_")]

    print(f"\n{sep}")
    print("  Multi-omics ML Integration — Summary")
    print(sep)
    print(f"\n  Master bins          : {len(master_df):>8,}")
    print(f"  Replicons            : {master_df['replicon'].nunique():>8,}")
    n_with_expr = master_df["expr_n_genes"].gt(0).sum()
    print(f"  Bins with expression : {n_with_expr:>8,}")
    print(f"  Methylation features : {len(methyl_cols):>8,}  (bin-level)")
    print(f"  Annotation features  : {len(annot_cols):>8,}  (gene-level)")

    print(f"\n  Model performance (test set):")
    print(f"  {sep2}")
    print(f"  {'Model':<20} {'RMSE':>10} {'R²':>10} {'Pearson r':>12}")
    print(f"  {sep2}")
    for m in [xgb_m, rf_m]:
        print(
            f"  {m['model']:<20} {m['rmse']:>10.4f} "
            f"{m['r_squared']:>10.4f} {m['pearson_r']:>12.4f}"
        )
    print(f"\n{sep}\n")


# ===========================================================================
# 8.  CLI
# ===========================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Integrate methylation, Hi-C contacts, and gene expression into a\n"
            "unified bin-level feature matrix, then train XGBoost + Random Forest\n"
            "models with SHAP interpretability."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  Basic:\n"
            "    python multiomics_ml.py \\\n"
            "        --annot-methyl annot_6mA.tsv \\\n"
            "        --bin-methyl   bin_6mA.bed \\\n"
            "        --deseq2       deseq2.tsv \\\n"
            "        --hic          contacts.bed \\\n"
            "        --output-dir   results/\n\n"
            "  Multiple methylation mod types:\n"
            "    python multiomics_ml.py \\\n"
            "        --annot-methyl annot_6mA.tsv annot_5mC.tsv \\\n"
            "        --bin-methyl   bin_6mA.bed bin_5mC.bed \\\n"
            "        --deseq2       deseq2.tsv \\\n"
            "        --hic          contacts.bed \\\n"
            "        --output-dir   results/ --resolution 10000\n"
        ),
    )

    req = p.add_argument_group("required inputs")
    req.add_argument(
        "--annot-methyl", nargs="+", required=True, metavar="FILE",
        help="Gene-level annotation methylation TSV (one or more files).",
    )
    req.add_argument(
        "--bin-methyl", nargs="+", required=True, metavar="FILE",
        help="Binned methylation BED (one or more files, e.g. per mod type).",
    )
    req.add_argument(
        "--deseq2", required=True, metavar="FILE",
        help="DESeq2 results TSV.",
    )
    req.add_argument(
        "--hic", required=True, metavar="FILE",
        help="Hi-C contact BED (chrom, start, end, contact_value).",
    )
    req.add_argument(
        "--output-dir", "-o", required=True, metavar="DIR",
        help="Directory for all outputs.",
    )

    opt = p.add_argument_group("optional parameters")
    opt.add_argument(
        "--resolution", type=int, default=None, metavar="INT",
        help=(
            "Override master bin size in bp. "
            "If omitted, the Hi-C bin size is used as the master resolution."
        ),
    )
    opt.add_argument(
        "--test-frac", type=float, default=0.2, metavar="FLOAT",
        help="Fraction of bins held out as test set [default: 0.2].",
    )
    opt.add_argument(
        "--seed", type=int, default=42, metavar="INT",
        help="Random seed [default: 42].",
    )
    opt.add_argument(
        "--format", default="pdf", choices=["pdf", "png", "svg"],
        metavar="EXT",
        help="Plot output format: pdf | png | svg [default: pdf].",
    )
    return p.parse_args()


# ===========================================================================
# 9.  Main
# ===========================================================================

def main() -> None:
    args = parse_args()

    out_dir = Path(args.output_dir).resolve()

    # ── Validate inputs ───────────────────────────────────────────────────────
    all_inputs = (
        [Path(f) for f in args.annot_methyl]
        + [Path(f) for f in args.bin_methyl]
        + [Path(args.deseq2), Path(args.hic)]
    )
    for f in all_inputs:
        if not f.is_file():
            logger.error("File not found: %s", f)
            sys.exit(1)

    if not HAS_XGB:
        logger.warning("xgboost not installed — only Random Forest will be trained.")
    if not HAS_SHAP:
        logger.warning("shap not installed — SHAP analysis will be skipped.")

    logger.info("=== Multi-omics ML Integration ===")

    # ── Load ──────────────────────────────────────────────────────────────────
    logger.info("Loading data...")
    annot_methyl = load_annot_methyl([Path(f) for f in args.annot_methyl])
    bin_methyl   = load_bin_methyl([Path(f) for f in args.bin_methyl])
    deseq2       = load_deseq2(Path(args.deseq2))
    hic_df       = load_hic(Path(args.hic))

    # ── Master bins ───────────────────────────────────────────────────────────
    logger.info("Building master bin grid...")
    master = build_master_bins(hic_df, args.resolution)

    # ── Aggregate each layer onto master bins ─────────────────────────────────
    logger.info("Aggregating data layers onto master bins...")
    bin_methyl_agg   = aggregate_bin_methyl(bin_methyl, master)
    hic_agg          = aggregate_hic(hic_df, master)
    expr_agg         = aggregate_expression(annot_methyl, deseq2, master)
    annot_methyl_agg = aggregate_annot_methyl(annot_methyl, master)

    # ── Build feature matrix ──────────────────────────────────────────────────
    logger.info("Building unified feature matrix...")
    master_df = build_feature_matrix(
        master, bin_methyl_agg, hic_agg, expr_agg, annot_methyl_agg
    )

    # Save the full matrix before subsetting to training data
    out_dir.mkdir(parents=True, exist_ok=True)
    master_df.to_csv(out_dir / "master_bins.tsv", sep="\t", index=False)
    logger.info("Saved: master_bins.tsv  (%d bins x %d columns)",
                len(master_df), len(master_df.columns))

    # ── Prepare ML data ───────────────────────────────────────────────────────
    logger.info("Preparing ML dataset...")
    X, y, feature_names = prepare_ml_data(master_df)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_frac, random_state=args.seed
    )
    logger.info("Train: %d  |  Test: %d", len(X_train), len(X_test))

    # ── Train XGBoost ─────────────────────────────────────────────────────────
    xgb_model, xgb_metrics = None, {"model": "XGBoost", "rmse": np.nan,
                                     "r_squared": np.nan, "pearson_r": np.nan}
    if HAS_XGB:
        logger.info("Training XGBoost...")
        xgb_model, xgb_metrics = train_xgboost(
            X_train, y_train, X_test, y_test, seed=args.seed
        )

    # ── Train Random Forest ───────────────────────────────────────────────────
    logger.info("Training Random Forest...")
    rf_model, rf_metrics = train_rf(
        X_train, y_train, X_test, y_test, seed=args.seed
    )

    # ── SHAP ─────────────────────────────────────────────────────────────────
    shap_vals = None
    if HAS_XGB and HAS_SHAP and xgb_model is not None:
        logger.info("Computing SHAP values (full training set)...")
        shap_vals = compute_shap(xgb_model, X_train)

    # ── Save ──────────────────────────────────────────────────────────────────
    logger.info("Saving outputs to %s ...", out_dir)
    save_outputs(
        out_dir, master_df, X_train, y_train,
        xgb_model, rf_model,
        xgb_metrics, rf_metrics,
        shap_vals, args.format,
    )

    # ── Summary ───────────────────────────────────────────────────────────────
    print_summary(master_df, xgb_metrics, rf_metrics)


if __name__ == "__main__":
    main()