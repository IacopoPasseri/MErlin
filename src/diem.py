#!/usr/bin/env python3
"""
 Integrates DESeq2 differential expression data with Oxford Nanopore-style
 methylation basecalling data to:
   1. Identify features that are both differentially expressed (DE) AND
      differentially methylated (DM) — the "dual-hit" set.
   2. Quantify Spearman & Pearson correlation between methylation level
      (norm_count) and log2FoldChange, with permutation-based significance.
   3. Build a Random Forest classifier to predict DE status from methylation
      features, with LOOCV (n<200) or stratified 5-fold CV (n>=200).
   4. Run univariate logistic regression to estimate the methylation odds
      ratio for DE, with BH-corrected q-values.
   5. Produce a multi-panel figure summarising all results.

 Input files:
   --geneexp    : DESeq2 output  CSV or TSV
                  (columns: deg, baseMean, log2FoldChange, lfcSE, stat,
                             pvalue, padj)
   --methylation: Basecalling methylation TSV
                  (columns: locus_tag, mod_label, norm_count, raw_count, …)

 Usage:
   python diem.py \
       --geneexp    geneexp_data.csv \
       --methylation methylation_data.tsv \
       --padj_thr   0.05 \
       --lfc_thr    1.0 \
       --meth_thr   0.004 \
       --outdir     results
=============================================================================
"""

import argparse
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

from scipy import stats
from scipy.stats import pearsonr, spearmanr

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (LeaveOneOut, StratifiedKFold,
                                     cross_val_predict, cross_val_score)
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (confusion_matrix, roc_auc_score, roc_curve,
                              ConfusionMatrixDisplay)

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# ── colour palette ─────────────────────────────────────────────────────────
PAL = {
    "DE":      "#e74c3c",
    "nonDE":   "#3498db",
    "DM":      "#f39c12",
    "dual":    "#9b59b6",
    "neutral": "#95a5a6",
}

# ═══════════════════════════════════════════════════════════════════════════
# 1.  DATA LOADING & VALIDATION
# ═══════════════════════════════════════════════════════════════════════════

def _read_table(path: str) -> pd.DataFrame:
    """Auto-detect CSV vs TSV from extension."""
    if path.lower().endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def load_and_validate(ge_path: str, me_path: str):
    ge = _read_table(ge_path)
    me = pd.read_csv(me_path, sep="\t")

    required_ge = {"deg", "log2FoldChange", "pvalue", "padj", "baseMean"}
    required_me = {"locus_tag", "norm_count", "raw_count", "mod_label"}
    missing_ge = required_ge - set(ge.columns)
    missing_me = required_me - set(me.columns)
    if missing_ge:
        raise ValueError(f"Gene expression table missing columns: {missing_ge}")
    if missing_me:
        raise ValueError(f"Methylation table missing columns: {missing_me}")

    ge = ge.dropna(subset=["padj", "log2FoldChange"])
    me = me.dropna(subset=["norm_count"])

    print(f"[INFO] Gene expression rows   : {len(ge):>6}  (unique loci: {ge['deg'].nunique()})")
    print(f"[INFO] Methylation rows       : {len(me):>6}  (unique loci: {me['locus_tag'].nunique()})")
    print(f"[INFO] Modification types     : {me['mod_label'].unique().tolist()}")

    overlap = set(ge["deg"]) & set(me["locus_tag"])
    print(f"[INFO] Overlapping loci       : {len(overlap)}")

    return ge, me


# ═══════════════════════════════════════════════════════════════════════════
# 2.  FEATURE ENGINEERING  –  pivot & aggregate methylation per locus
# ═══════════════════════════════════════════════════════════════════════════

def aggregate_methylation(me: pd.DataFrame) -> pd.DataFrame:
    """
    For each locus_tag × mod_label produce:
      total_norm_count – sum of normalised methylation counts across sites
      mean_norm_count  – mean normalised count
      raw_count_total  – total raw methylation events
      n_sites          – number of independent sites detected
    Then pivot so each mod_label becomes a column prefix (one row per locus).
    """
    agg = (
        me.groupby(["locus_tag", "mod_label"])
        .agg(
            total_norm_count=("norm_count", "sum"),
            mean_norm_count =("norm_count", "mean"),
            raw_count_total =("raw_count",  "sum"),
            n_sites         =("norm_count", "count"),
        )
        .reset_index()
    )

    pivot = agg.pivot(index="locus_tag", columns="mod_label")
    pivot.columns = ["_".join(c).strip() for c in pivot.columns]
    pivot = pivot.reset_index().fillna(0)

    print(f"\n[METHYLATION PIVOT]  shape: {pivot.shape}")
    print(f"  Columns: {pivot.columns.tolist()}")
    return pivot


# ═══════════════════════════════════════════════════════════════════════════
# 3.  DEFINE DE AND DM GROUPS
# ═══════════════════════════════════════════════════════════════════════════

def define_groups(ge, me_pivot, padj_thr, lfc_thr, meth_thr):
    """
    Merge on locus_tag (deg ↔ locus_tag) and annotate:
      DE   : padj < padj_thr  AND  |log2FoldChange| >= lfc_thr
      DM   : max total_norm_count across mod types > meth_thr
      dual : DE ∩ DM
    Duplicate deg entries are collapsed by keeping the minimum-padj row.
    """
    ge_agg = ge.sort_values("padj").groupby("deg", as_index=False).first()
    ge_agg = ge_agg.rename(columns={"deg": "locus_tag"})

    merged = ge_agg.merge(me_pivot, on="locus_tag", how="inner")

    merged["is_DE"] = (
        (merged["padj"] < padj_thr) &
        (merged["log2FoldChange"].abs() >= lfc_thr)
    ).astype(int)

    meth_total_cols = [c for c in merged.columns if c.startswith("total_norm_count")]
    merged["max_meth"] = merged[meth_total_cols].max(axis=1) if meth_total_cols else 0

    merged["is_DM"]    = (merged["max_meth"] > meth_thr).astype(int)
    merged["dual_hit"] = ((merged["is_DE"] == 1) & (merged["is_DM"] == 1)).astype(int)

    n_de   = int(merged["is_DE"].sum())
    n_dm   = int(merged["is_DM"].sum())
    n_dual = int(merged["dual_hit"].sum())

    print(f"\n[GROUPS]  n = {len(merged):,} overlapping loci")
    print(f"  DE  (padj<{padj_thr}, |LFC|≥{lfc_thr})      : {n_de:,}")
    print(f"  DM  (max_norm_count > {meth_thr})             : {n_dm:,}")
    print(f"  Dual-hit (DE ∩ DM)                            : {n_dual:,}")

    return merged


# ═══════════════════════════════════════════════════════════════════════════
# 4.  CORRELATION ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

def correlation_analysis(df: pd.DataFrame) -> dict:
    """Pearson & Spearman between max methylation level and log2FoldChange."""
    x = df["max_meth"].values
    y = df["log2FoldChange"].values

    r_p, pval_p = pearsonr(x, y)
    r_s, pval_s = spearmanr(x, y)

    # Permutation p-value for Spearman (5 000 shuffles)
    n_perm  = 5000
    null_rs = np.array([spearmanr(np.random.permutation(x), y)[0]
                        for _ in range(n_perm)])
    perm_pval = float(np.mean(np.abs(null_rs) >= abs(r_s)))

    # Per-mod Spearman breakdown
    mod_corrs = {}
    for col in [c for c in df.columns if c.startswith("total_norm_count_")]:
        mod = col.replace("total_norm_count_", "")
        r_m, p_m = spearmanr(df[col].values, y)
        mod_corrs[mod] = (r_m, p_m)

    results = {
        "pearson_r": r_p, "pearson_pval": pval_p,
        "spearman_r": r_s, "spearman_pval": pval_s,
        "spearman_perm_pval": perm_pval,
        "mod_corrs": mod_corrs,
    }

    print("\n[CORRELATION]  methylation (max_norm_count) vs log2FoldChange")
    print(f"  Pearson  r = {r_p:+.4f}  p = {pval_p:.4e}")
    print(f"  Spearman ρ = {r_s:+.4f}  p = {pval_s:.4e}  (perm p = {perm_pval:.4f})")
    for mod, (r_m, p_m) in mod_corrs.items():
        print(f"  Spearman [{mod}] ρ = {r_m:+.4f}  p = {p_m:.4e}")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# 5.  LOGISTIC REGRESSION  (univariate per feature)
# ═══════════════════════════════════════════════════════════════════════════

def logistic_regression_analysis(df: pd.DataFrame, feature_cols: list) -> pd.DataFrame:
    """
    Univariate logistic regression: predict is_DE from each methylation feature.
    Returns a table with coefficient, OR, 95 % CI, p-value, BH q-value.
    """
    records = []
    for col in feature_cols:
        vals = df[col].astype(float)
        if vals.std() == 0:
            continue
        X = sm.add_constant(vals)
        y = df["is_DE"]
        if y.nunique() < 2:
            continue
        try:
            model = sm.Logit(y, X).fit(disp=0, maxiter=200)
            coef  = model.params[col]
            ci    = model.conf_int().loc[col]
            pval  = model.pvalues[col]
            records.append({
                "feature": col,
                "coef":    coef,
                "OR":      np.exp(coef),
                "CI_low":  np.exp(ci[0]),
                "CI_high": np.exp(ci[1]),
                "pvalue":  pval,
            })
        except Exception as e:
            print(f"  [WARN] LR failed for {col}: {e}")

    if not records:
        print("[WARN] No logistic regression models converged.")
        return pd.DataFrame()

    res = pd.DataFrame(records)
    _, res["qvalue"], _, _ = multipletests(res["pvalue"], method="fdr_bh")
    res = res.sort_values("pvalue").reset_index(drop=True)

    print("\n[LOGISTIC REGRESSION]  univariate results")
    print(res.to_string(index=False, float_format="{:.4f}".format))
    return res


# ═══════════════════════════════════════════════════════════════════════════
# 6.  RANDOM FOREST  (LOOCV if n<200, else stratified 5-fold)
# ═══════════════════════════════════════════════════════════════════════════

def random_forest_cv(df: pd.DataFrame, feature_cols: list) -> dict:
    """
    Train a Random Forest to predict is_DE from methylation features.

    Model choice rationale
    ─────────────────────
    • Random Forest handles non-linear methylation–expression relationships
      without distributional assumptions.
    • Robust to multi-collinearity between 6mA and 5mC features.
    • class_weight='balanced' corrects for DE/non-DE imbalance.
    • LOOCV for n<200 (maximises training data); stratified 5-fold for
      larger datasets (computationally tractable, still unbiased).
    """
    X = df[feature_cols].fillna(0).astype(float).values
    y = df["is_DE"].values

    if y.sum() < 2 or (y == 0).sum() < 2:
        print("[WARN] Insufficient class diversity for Random Forest.")
        return {}

    scaler = StandardScaler()
    X_sc   = scaler.fit_transform(X)

    rf = RandomForestClassifier(
        n_estimators=500,
        max_depth=None,
        min_samples_leaf=2,
        class_weight="balanced",
        random_state=42,
        n_jobs=-1,
    )

    n = len(df)
    if n < 200:
        cv = LeaveOneOut()
        cv_name = "LOOCV"
    else:
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        cv_name = "5-fold stratified CV"

    print(f"\n[RANDOM FOREST]  CV strategy: {cv_name}  (n={n:,})")

    y_pred_proba = cross_val_predict(rf, X_sc, y, cv=cv, method="predict_proba")[:, 1]
    y_pred       = (y_pred_proba >= 0.5).astype(int)

    accuracy = float(np.mean(y_pred == y))
    try:
        auc = float(roc_auc_score(y, y_pred_proba))
    except Exception:
        auc = float("nan")

    cm = confusion_matrix(y, y_pred)

    # Fit on full data for feature importances
    rf.fit(X_sc, y)
    importances = pd.Series(rf.feature_importances_,
                            index=feature_cols).sort_values(ascending=False)

    print(f"  Accuracy : {accuracy:.4f}")
    print(f"  AUC-ROC  : {auc:.4f}")
    print(f"  Confusion matrix:\n{cm}")
    print(f"  Feature importances:\n{importances.to_string()}")

    return {
        "model": rf, "scaler": scaler, "cv_name": cv_name,
        "accuracy": accuracy, "auc": auc,
        "y_true": y, "y_pred": y_pred, "y_pred_proba": y_pred_proba,
        "confusion_matrix": cm, "feature_importances": importances,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 7.  FIGURE  (9-panel)
# ═══════════════════════════════════════════════════════════════════════════

def _short(col: str) -> str:
    """Shorten a feature column name for axis labels."""
    return (col.replace("total_norm_count_", "Σnorm_")
               .replace("mean_norm_count_",  "μnorm_")
               .replace("raw_count_total_",  "raw_")
               .replace("n_sites_",          "sites_"))


def make_figure(df, corr, rf_res, lr_res, outdir):
    fig = plt.figure(figsize=(20, 15), facecolor="#fafafa")
    fig.suptitle("Methylation – Expression Integrative Analysis",
                 fontsize=16, fontweight="bold", y=0.99, color="#2c3e50")

    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

    # colour per point
    def point_color(row):
        if row["dual_hit"]:   return PAL["dual"]
        if row["is_DE"]:      return PAL["DE"]
        if row["is_DM"]:      return PAL["DM"]
        return PAL["neutral"]
    cmap_vec = df.apply(point_color, axis=1)

    # ── A: Volcano ─────────────────────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    neg_log_p = -np.log10(df["padj"].clip(lower=1e-300))
    ax_a.scatter(df["log2FoldChange"], neg_log_p,
                 c=cmap_vec, s=8, alpha=0.6, linewidths=0)
    ax_a.axvline(0,  lw=0.8, ls="--", color="#555")
    ax_a.axhline(-np.log10(0.05), lw=0.8, ls=":", color="#555")
    ax_a.set_xlabel("log₂ Fold Change", fontsize=9)
    ax_a.set_ylabel("−log₁₀(padj)", fontsize=9)
    ax_a.set_title("A  Volcano Plot", fontsize=10, fontweight="bold")
    _legend_patches(ax_a)

    # ── B: Methylation vs LFC ──────────────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.scatter(df["max_meth"], df["log2FoldChange"],
                 c=cmap_vec, s=8, alpha=0.5, linewidths=0)
    m, b_val, *_ = stats.linregress(df["max_meth"], df["log2FoldChange"])
    xline = np.linspace(df["max_meth"].min(), df["max_meth"].max(), 200)
    ax_b.plot(xline, m * xline + b_val, lw=1.8, color="#2c3e50", ls="--")
    r_lab = (f"Pearson r={corr['pearson_r']:+.3f} (p={corr['pearson_pval']:.2e})\n"
             f"Spearman ρ={corr['spearman_r']:+.3f} (p={corr['spearman_pval']:.2e})")
    ax_b.set_xlabel("Max normalised methylation count", fontsize=9)
    ax_b.set_ylabel("log₂ Fold Change", fontsize=9)
    ax_b.set_title(f"B  Methylation vs LFC\n{r_lab}", fontsize=9, fontweight="bold")

    # ── C: Group summary bar ──────────────────────────────────────────────
    ax_c = fig.add_subplot(gs[0, 2])
    cats   = ["Overlap", "DE", "DM", "Dual-hit"]
    counts = [len(df), int(df["is_DE"].sum()), int(df["is_DM"].sum()), int(df["dual_hit"].sum())]
    colors = [PAL["neutral"], PAL["DE"], PAL["DM"], PAL["dual"]]
    bars   = ax_c.bar(cats, counts, color=colors, edgecolor="white", linewidth=0.8)
    for b_, v in zip(bars, counts):
        ax_c.text(b_.get_x() + b_.get_width() / 2, v + max(counts)*0.01,
                  f"{v:,}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax_c.set_ylabel("Number of loci", fontsize=9)
    ax_c.set_title("C  Feature summary", fontsize=10, fontweight="bold")
    ax_c.set_ylim(0, max(counts) * 1.18)
    ax_c.tick_params(axis="x", labelsize=8)

    # ── D: Boxplot methylation by mod type & DE status ────────────────────
    ax_d = fig.add_subplot(gs[1, 0])
    meth_total_cols = [c for c in df.columns if c.startswith("total_norm_count")]
    if len(meth_total_cols) >= 1:
        melt = df[["is_DE"] + meth_total_cols].melt(
            id_vars="is_DE", var_name="mod", value_name="norm_count")
        melt["mod"] = melt["mod"].str.replace("total_norm_count_", "", regex=False)
        melt["DE Status"] = melt["is_DE"].map({1: "DE", 0: "non-DE"})
        sns.boxplot(data=melt, x="mod", y="norm_count", hue="DE Status",
                    palette={"DE": PAL["DE"], "non-DE": PAL["nonDE"]},
                    ax=ax_d, fliersize=2, linewidth=0.8, showfliers=False)
        # Mann-Whitney U p-values per mod
        for i, mod in enumerate(melt["mod"].unique()):
            sub = melt[melt["mod"] == mod]
            grp0 = sub.loc[sub["DE Status"] == "non-DE", "norm_count"]
            grp1 = sub.loc[sub["DE Status"] == "DE",     "norm_count"]
            if len(grp0) > 1 and len(grp1) > 1:
                _, p = stats.mannwhitneyu(grp0, grp1, alternative="two-sided")
                y_max = sub["norm_count"].quantile(0.97)
                ax_d.text(i, y_max * 1.02, f"p={p:.2e}", ha="center",
                          fontsize=7, color="#333")
        ax_d.set_xlabel("Modification type", fontsize=9)
        ax_d.set_ylabel("Total norm. count per locus", fontsize=9)
        ax_d.set_title("D  Methylation by mod & DE status", fontsize=10, fontweight="bold")
        ax_d.legend(fontsize=7)

    # ── E: Feature importances ────────────────────────────────────────────
    ax_e = fig.add_subplot(gs[1, 1])
    if rf_res and "feature_importances" in rf_res:
        fi = rf_res["feature_importances"]
        labels = [_short(i) for i in fi.index]
        colors_fi = [PAL["dual"] if v == fi.max() else PAL["nonDE"] for v in fi.values]
        ax_e.barh(labels[::-1], fi.values[::-1], color=colors_fi[::-1], edgecolor="white")
        ax_e.set_xlabel("Mean Decrease Impurity", fontsize=9)
        cv_tag = rf_res.get("cv_name", "CV")
        ax_e.set_title(
            f"E  RF Feature Importances\n({cv_tag}  Acc={rf_res['accuracy']:.3f}  AUC={rf_res['auc']:.3f})",
            fontsize=9, fontweight="bold")
    else:
        ax_e.text(0.5, 0.5, "RF not fitted", ha="center", va="center",
                  transform=ax_e.transAxes, color="grey")
        ax_e.set_title("E  RF Feature Importances", fontsize=10, fontweight="bold")

    # ── F: Confusion matrix ───────────────────────────────────────────────
    ax_f = fig.add_subplot(gs[1, 2])
    if rf_res and "confusion_matrix" in rf_res:
        cv_tag = rf_res.get("cv_name", "CV")
        disp = ConfusionMatrixDisplay(rf_res["confusion_matrix"],
                                      display_labels=["non-DE", "DE"])
        disp.plot(ax=ax_f, colorbar=False, cmap="Blues")
        ax_f.set_title(f"F  Confusion Matrix ({cv_tag})", fontsize=10, fontweight="bold")
    else:
        ax_f.text(0.5, 0.5, "RF not fitted", ha="center", va="center",
                  transform=ax_f.transAxes, color="grey")
        ax_f.set_title("F  Confusion Matrix", fontsize=10, fontweight="bold")

    # ── G: ROC curve ──────────────────────────────────────────────────────
    ax_g = fig.add_subplot(gs[2, 0])
    if rf_res and "y_pred_proba" in rf_res:
        fpr, tpr, _ = roc_curve(rf_res["y_true"], rf_res["y_pred_proba"])
        ax_g.plot(fpr, tpr, lw=2, color=PAL["dual"],
                  label=f"AUC = {rf_res['auc']:.3f}")
        ax_g.plot([0, 1], [0, 1], lw=1, ls="--", color="grey")
        ax_g.set_xlabel("False Positive Rate", fontsize=9)
        ax_g.set_ylabel("True Positive Rate", fontsize=9)
        cv_tag = rf_res.get("cv_name", "CV")
        ax_g.set_title(f"G  ROC Curve ({cv_tag})", fontsize=10, fontweight="bold")
        ax_g.legend(fontsize=9)
    else:
        ax_g.text(0.5, 0.5, "RF not fitted", ha="center", va="center",
                  transform=ax_g.transAxes, color="grey")
        ax_g.set_title("G  ROC Curve", fontsize=10, fontweight="bold")

    # ── H: LR forest plot ─────────────────────────────────────────────────
    ax_h = fig.add_subplot(gs[2, 1])
    if lr_res is not None and not lr_res.empty:
        labels_h = [_short(f) for f in lr_res["feature"]]
        y_pos = range(len(lr_res))
        # clip extreme CI for display
        or_  = lr_res["OR"].clip(1e-4, 1e4)
        cil  = (lr_res["OR"] - lr_res["CI_low"]).clip(0)
        cih  = (lr_res["CI_high"] - lr_res["OR"]).clip(0)
        ax_h.errorbar(or_, list(y_pos),
                      xerr=[cil, cih],
                      fmt="o", color=PAL["DE"], ecolor="#aaa",
                      capsize=4, ms=6, linewidth=1.2)
        ax_h.axvline(1.0, lw=1, ls="--", color="#555")
        ax_h.set_xscale("log")
        ax_h.set_yticks(list(y_pos))
        ax_h.set_yticklabels(labels_h, fontsize=8)
        ax_h.set_xlabel("Odds Ratio (95 % CI, log scale)", fontsize=9)
        ax_h.set_title("H  Logistic Regression\nOR forest plot", fontsize=9, fontweight="bold")
        # mark significant features
        for i, row in lr_res.iterrows():
            if row["qvalue"] < 0.05:
                ax_h.text(or_.iloc[i] * 1.05, i, "*", fontsize=11, color=PAL["DE"], va="center")
    else:
        ax_h.text(0.5, 0.5, "LR not fitted", ha="center", va="center",
                  transform=ax_h.transAxes, color="grey")
        ax_h.set_title("H  Logistic Regression", fontsize=10, fontweight="bold")

    # ── I: Methylation heatmap (top 40 dual-hit loci) ─────────────────────
    ax_i = fig.add_subplot(gs[2, 2])
    heat_cols = [c for c in df.columns if c.startswith("total_norm_count") or
                 c.startswith("mean_norm_count")]
    if heat_cols:
        plot_df = df[df["dual_hit"] == 1] if df["dual_hit"].sum() >= 5 else df
        plot_df = plot_df.set_index("locus_tag")[heat_cols]
        # Limit to top 40 by max methylation
        if len(plot_df) > 40:
            plot_df = plot_df.loc[plot_df.max(axis=1).nlargest(40).index]
        plot_df.columns = [_short(c) for c in plot_df.columns]
        annot = len(plot_df) <= 40
        sns.heatmap(plot_df, ax=ax_i, cmap="YlOrRd", linewidths=0.2,
                    linecolor="white", annot=annot,
                    fmt=".4f" if annot else "", annot_kws={"size": 6},
                    cbar_kws={"shrink": 0.7})
        title_tag = "dual-hit" if df["dual_hit"].sum() >= 5 else "all overlapping"
        ax_i.set_title(f"I  Methylation Heatmap\n(top {len(plot_df)} {title_tag} loci)",
                       fontsize=9, fontweight="bold")
        ax_i.tick_params(axis="x", rotation=30, labelsize=8)
        ax_i.tick_params(axis="y", rotation=0,  labelsize=6)
    else:
        ax_i.text(0.5, 0.5, "No methylation columns", ha="center", va="center",
                  transform=ax_i.transAxes, color="grey")
        ax_i.set_title("I  Methylation Heatmap", fontsize=10, fontweight="bold")

    out_fig = os.path.join(outdir, "meth_expr_integration.png")
    fig.savefig(out_fig, dpi=180, bbox_inches="tight")
    print(f"\n[FIGURE] Saved → {out_fig}")
    plt.close(fig)


def _legend_patches(ax):
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor=PAL["DE"],      label="DE only"),
        Patch(facecolor=PAL["DM"],      label="DM only"),
        Patch(facecolor=PAL["dual"],    label="Dual-hit"),
        Patch(facecolor=PAL["neutral"], label="Neither"),
    ], fontsize=6, loc="upper right")


# ═══════════════════════════════════════════════════════════════════════════
# 8.  SUMMARY TABLES
# ═══════════════════════════════════════════════════════════════════════════

def save_summary(df, corr, rf_res, lr_res, outdir):
    df.to_csv(os.path.join(outdir, "integrated_table.tsv"), sep="\t", index=False)
    df[df["dual_hit"] == 1].to_csv(
        os.path.join(outdir, "dual_hit_features.tsv"), sep="\t", index=False)

    lines = [
        "=== Correlation  (methylation max_norm_count vs log2FoldChange) ===",
        f"Pearson  r={corr['pearson_r']:+.4f}  p={corr['pearson_pval']:.4e}",
        f"Spearman ρ={corr['spearman_r']:+.4f}  p={corr['spearman_pval']:.4e}"
        f"  (perm p={corr['spearman_perm_pval']:.4f})",
        "",
        "=== Per-modification Spearman correlations ===",
    ]
    for mod, (r_m, p_m) in corr.get("mod_corrs", {}).items():
        lines.append(f"  [{mod}]  ρ={r_m:+.4f}  p={p_m:.4e}")

    lines += ["", "=== Random Forest ==="]
    if rf_res:
        lines += [
            f"CV strategy : {rf_res.get('cv_name','CV')}",
            f"Accuracy    : {rf_res['accuracy']:.4f}",
            f"AUC-ROC     : {rf_res['auc']:.4f}",
        ]
    else:
        lines.append("Not fitted.")

    with open(os.path.join(outdir, "statistics_summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    if lr_res is not None and not lr_res.empty:
        lr_res.to_csv(os.path.join(outdir, "logistic_regression_results.tsv"),
                      sep="\t", index=False)

    print(f"[OUTPUT] Tables written to: {outdir}/")


# ═══════════════════════════════════════════════════════════════════════════
# 9.  MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="DIEM: Differential Integration of Expression and Methylation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--geneexp",     required=True,
                        help="DESeq2 CSV or TSV file")
    parser.add_argument("--methylation", required=True,
                        help="Methylation basecalling TSV")
    parser.add_argument("--padj_thr",    type=float, default=0.05,
                        help="Adjusted p-value cut-off for DE")
    parser.add_argument("--lfc_thr",     type=float, default=1.0,
                        help="|log2FoldChange| cut-off for DE")
    parser.add_argument("--meth_thr",    type=float, default=0.004,
                        help="Max normalised count cut-off for DM (≈75th pct of real data)")
    parser.add_argument("--outdir",      default="results",
                        help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("=" * 60)
    print("  DIEM: Differential Integration of Expression and Methylation")
    print("=" * 60)

    ge, me     = load_and_validate(args.geneexp, args.methylation)
    me_pivot   = aggregate_methylation(me)
    df         = define_groups(ge, me_pivot, args.padj_thr, args.lfc_thr, args.meth_thr)

    corr       = correlation_analysis(df)

    feature_cols = [c for c in df.columns if any(
        c.startswith(p) for p in
        ["total_norm_count", "mean_norm_count", "raw_count_total", "n_sites"]
    )]
    print(f"\n[ML] Feature columns ({len(feature_cols)}): {feature_cols}")

    lr_res     = logistic_regression_analysis(df, feature_cols)
    rf_res     = random_forest_cv(df, feature_cols)

    make_figure(df, corr, rf_res, lr_res, args.outdir)
    save_summary(df, corr, rf_res, lr_res, args.outdir)

    print("\n[DONE] Pipeline complete.\n")


if __name__ == "__main__":
    main()
