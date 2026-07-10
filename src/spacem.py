#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spatial_oriter_distribution.py

Describe (statistics) and visualise (plots) the spatial distribution of DNA
methylation calls and genomic features relative to the replication origin (oriC)
and terminus (ter), to assess whether their genomic arrangement is "spatially
functional" -- i.e. organised by the replication geometry rather than randomly.

Takes the SAME inputs as the cross-referencing tool:
  --gff (GFF3), --fasta (reference), --methylation (bedMethyl or GFF3),
  and optionally --motifs (to add a per-motif spatial track).

WHAT IT COMPUTES
----------------
1. oriC / ter per replicon from cumulative GC-skew (origin = global minimum of
   the cumulative (G-C) curve, terminus = global maximum). Both can be overridden
   (--oric / --ter). Each replicon is treated as circular.
2. An ori->ter coordinate for every position: replication runs bidirectionally
   from oriC to ter along two replichores; each base is assigned to a replichore
   and a normalised distance from oriC (0 = oriC, 1 = ter). A signed coordinate
   in [-1, +1] keeps the two replichores distinguishable.
3. Gene orientation relative to the replication fork (co-oriented = on the
   leading strand). Strong co-orientation bias is a classic replication-associated
   (spatially functional) signal.
4. Statistics for BOTH layers:
     * gene co-orientation fraction + binomial test vs 0.5,
     * spatial non-uniformity of methylation and of genes along the ori->ter axis
       (KS test of the folded distance-from-oriC against Uniform(0,1); mean
       distance < 0.5 => ori-biased, > 0.5 => ter-biased),
     * methylation-LEVEL gradient: Spearman correlation of per-call fraction
       modified vs distance from oriC (tests replication-timing-like gradients,
       e.g. CcrM/GANTC in alphaproteobacteria).

OUTPUTS (--outdir, prefixed with --prefix)
------------------------------------------
  <prefix>.oriter_stats.tsv     one row per replicon (+ a genome row) with all
                                statistics above.
  <prefix>.spatial_bins.tsv     per-replicon binned profiles (gene density,
                                methylation density, mean fraction, co-orientation).
  <prefix>.<replicon>.spatial.<fmt>      4-panel figure per replicon (with --plots)
  <prefix>.<replicon>.motif_spatial.<fmt>  per-motif methylation track (if --motifs)

Interpretation note: methylation calls only exist where a base was assayed with
coverage, so methylation DENSITY partly reflects assay/sequence context; the
LEVEL gradient (fraction vs distance) and gene co-orientation are the more robust
"spatially functional" readouts. Association is not causation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

__version__ = "1.1.0"


# --------------------------------------------------------------------------- #
# Small IO + parsing helpers (self-contained)
# --------------------------------------------------------------------------- #
def smart_open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


_VERSION_RE = re.compile(r"\.\d+$")


def strip_version(s: str) -> str:
    return _VERSION_RE.sub("", s)


@dataclass
class Feature:
    seqid: str
    start: int   # 0-based inclusive
    end: int     # 0-based inclusive
    strand: str
    ftype: str
    fid: str


@dataclass
class MethCall:
    seqid: str
    pos: int     # 0-based
    strand: str
    mod_code: str
    coverage: Optional[float]
    frac: Optional[float]


def parse_gff3(path, feature_types, id_attr) -> Tuple[List[Feature], Dict[str, str]]:
    features: List[Feature] = []
    embedded: Dict[str, str] = {}
    in_fasta = False
    cid = None
    buf: List[str] = []
    auto = 0
    with smart_open(path) as fh:
        for line in fh:
            if in_fasta:
                if line.startswith(">"):
                    if cid is not None:
                        embedded[cid] = "".join(buf).upper()
                    cid = line[1:].strip().split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
                continue
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            seqid, _s, ftype, start, end, _sc, strand, _ph, attr = c[:9]
            if feature_types and ftype not in feature_types:
                continue
            try:
                s0, e0 = int(start) - 1, int(end) - 1
            except ValueError:
                continue
            attrs = {}
            for kv in attr.strip().split(";"):
                kv = kv.strip()
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = v.strip()
            fid = attrs.get(id_attr) or attrs.get("locus_tag") or attrs.get("Name")
            if not fid:
                auto += 1
                fid = f"{ftype}_{seqid}_{s0+1}_{auto}"
            features.append(Feature(seqid, s0, e0, strand, ftype, fid))
    if in_fasta and cid is not None:
        embedded[cid] = "".join(buf).upper()
    return features, embedded


def load_fasta(path) -> Dict[str, str]:
    g: Dict[str, str] = {}
    cid = None
    buf: List[str] = []
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cid is not None:
                    g[cid] = "".join(buf).upper()
                cid = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cid is not None:
        g[cid] = "".join(buf).upper()
    return g


def _norm_frac(v: float) -> float:
    return v / 100.0 if v > 1.0 else v


def parse_methylation(path, fmt, keep_codes) -> List[MethCall]:
    low = path.lower()
    if fmt == "auto":
        if low.endswith((".bed", ".bed.gz", ".bedmethyl", ".bedmethyl.gz")):
            fmt = "bedmethyl"
        elif low.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
            fmt = "gff"
        else:
            raise ValueError("cannot infer methylation format; use --meth-format")
    calls: List[MethCall] = []
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            if fmt == "bedmethyl":
                c = line.rstrip("\n").split("\t")
                if len(c) < 6:
                    c = line.split()
                if len(c) < 6:
                    continue
                try:
                    pos = int(c[1])
                except ValueError:
                    continue
                code = c[3] if len(c) > 3 else "."
                if keep_codes and code not in keep_codes:
                    continue
                strand = c[5] if c[5] in ("+", "-", ".") else "."
                cov = frac = None
                if len(c) >= 11:
                    try:
                        cov = float(c[9])
                    except ValueError:
                        pass
                    try:
                        frac = _norm_frac(float(c[10]))
                    except ValueError:
                        pass
                calls.append(MethCall(c[0], pos, strand, code, cov, frac))
            else:  # gff
                c = line.rstrip("\n").split("\t")
                if len(c) < 8:
                    continue
                code = c[2]
                if keep_codes and code not in keep_codes:
                    continue
                try:
                    pos = int(c[3]) - 1
                except ValueError:
                    continue
                strand = c[6] if c[6] in ("+", "-", ".") else "."
                attrs = {}
                if len(c) >= 9:
                    for kv in c[8].split(";"):
                        if "=" in kv:
                            k, v = kv.split("=", 1)
                            attrs[k.strip()] = v.strip()
                cov = frac = None
                for k in ("coverage", "cov", "Nvalid_cov"):
                    if k in attrs:
                        try:
                            cov = float(attrs[k]); break
                        except ValueError:
                            pass
                for k in ("frac", "frac_modified", "fraction", "percent_modified"):
                    if k in attrs:
                        try:
                            frac = _norm_frac(float(attrs[k])); break
                        except ValueError:
                            pass
                calls.append(MethCall(c[0], pos, strand, code, cov, frac))
    return calls


# --------------------------------------------------------------------------- #
# IUPAC + motif search (only used when --motifs is given)
# --------------------------------------------------------------------------- #
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
         "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]",
         "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
         "N": "[ACGT]"}
COMP = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G", "R": "Y", "Y": "R",
        "S": "S", "W": "W", "K": "M", "M": "K", "B": "V", "V": "B", "D": "H",
        "H": "D", "N": "N"}


def revcomp(m: str) -> str:
    return "".join(COMP[b] for b in reversed(m.upper()))


def motif_regex(m: str):
    return re.compile("(?=(" + "".join(IUPAC[b] for b in m.upper()) + "))")


def parse_motifs(path) -> List[str]:
    out = []
    with smart_open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                out.append(s.split()[0].upper())
    return out


def motif_covered_positions(genome, motifs) -> Dict[str, set]:
    """Per motif: set of (seqid, pos) covered by any occurrence (both strands)."""
    _occ, cov = find_motif_occurrences(genome, motifs)
    return cov


def find_motif_occurrences(genome, motifs):
    """Per motif, scan both strands of every replicon.

    Returns (occ, cov):
      occ[motif][seqid] -> list of occurrence start positions (0-based)
      cov[motif]        -> set of (seqid, pos) covered by any occurrence
    """
    occ: Dict[str, Dict[str, list]] = {m: defaultdict(list) for m in motifs}
    cov: Dict[str, set] = {m: set() for m in motifs}
    for m in motifs:
        L = len(m)
        fwd = motif_regex(m)
        rc = revcomp(m)
        rev = fwd if rc == m else motif_regex(rc)
        for seqid, seq in genome.items():
            cset = cov[m]
            olist = occ[m][seqid]
            for mt in fwd.finditer(seq):
                p = mt.start()
                olist.append(p)
                for i in range(p, p + L):
                    cset.add((seqid, i))
            for mt in rev.finditer(seq):
                p = mt.start()
                olist.append(p)
                for i in range(p, p + L):
                    cset.add((seqid, i))
    return occ, cov


# --------------------------------------------------------------------------- #
# oriC / ter from cumulative GC skew
# --------------------------------------------------------------------------- #
@dataclass
class SkewProfile:
    centers: List[int]
    cum: List[int]
    oric: int
    ter: int
    amplitude: int


def gc_skew(seq: str, window: int) -> SkewProfile:
    centers, deltas, cumlist = [], [], []
    cum = 0
    L = len(seq)
    for i in range(0, L, window):
        w = seq[i:i + window]
        d = w.count("G") - w.count("C")
        deltas.append(d)
        cum += d
        centers.append(i + len(w) // 2)
        cumlist.append(cum)
    n = len(deltas)
    # oriC = global minimum of the (linear) cumulative skew -- robust to the
    # arbitrary contig start.
    oric_idx = min(range(n), key=lambda k: cumlist[k])
    # ter = maximum of the cumulative skew RE-STARTED at oriC, walking the circle.
    # Using the linear argmax instead would lock onto the contig start, not the
    # true antipodal terminus.
    s = 0
    rot = []
    for j in range(n):
        s += deltas[(oric_idx + j) % n]
        rot.append(s)
    ter_idx = (oric_idx + max(range(n), key=lambda k: rot[k])) % n
    amplitude = max(rot) - min(rot)
    return SkewProfile(centers, cumlist, centers[oric_idx], centers[ter_idx],
                       amplitude)


# --------------------------------------------------------------------------- #
# ori->ter coordinate
# --------------------------------------------------------------------------- #
def ori_ter(pos: int, oric: int, ter: int, L: int):
    """Return (signed_norm in [-1,1], abs_norm in [0,1], arm 'R+'/'R-')."""
    d = (pos - oric) % L
    T = (ter - oric) % L
    if T == 0:
        T = L // 2
    if d <= T:
        x = d / T if T else 0.0
        return x, x, "R+"
    back = L - d
    denom = L - T
    x = back / denom if denom else 0.0
    return -x, x, "R-"


def co_oriented(strand: str, arm: str) -> Optional[bool]:
    """Leading strand is '+' on the R+ arm and '-' on the R- arm."""
    if strand not in ("+", "-"):
        return None
    lead = "+" if arm == "R+" else "-"
    return strand == lead


# --------------------------------------------------------------------------- #
# Stats helpers
# --------------------------------------------------------------------------- #
def try_scipy():
    try:
        import scipy.stats as st
        return st
    except ImportError:
        return None


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def fold_enrichment(signed_vals, nbins):
    """Histogram over [-1,1]; normalise each replichore half by its own mean so a
    spatially uniform distribution reads ~1.0. Returns (centers, fold)."""
    import numpy as np
    counts, edges = np.histogram(signed_vals, bins=nbins, range=(-1, 1))
    centers = (edges[:-1] + edges[1:]) / 2
    fold = counts.astype(float)
    neg = centers < 0
    pos = ~neg
    for mask in (neg, pos):
        mu = fold[mask].mean() if fold[mask].size and fold[mask].mean() else 1.0
        fold[mask] = fold[mask] / mu if mu else 0.0
    return centers, fold


def make_replicon_figure(outpath, repname, length, skew: SkewProfile,
                         g_signed, m_signed, m_absfrac, co_abs, nbins, plot_fmt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"{repname}  ({length/1e6:.2f} Mb)   spatial organisation vs oriC/ter",
                 fontsize=14, fontweight="bold")

    # (a) cumulative GC skew
    ax = axes[0, 0]
    xs = np.array(skew.centers) / 1e6
    ax.plot(xs, skew.cum, color="#34495e", lw=1.2)
    ax.axvline(skew.oric / 1e6, color="#27ae60", lw=1.6, label="oriC")
    ax.axvline(skew.ter / 1e6, color="#c0392b", lw=1.6, label="ter")
    ax.set_xlabel("genomic position (Mb)")
    ax.set_ylabel("cumulative (G$-$C)")
    ax.set_title("(a) Cumulative GC skew -> oriC/ter")
    ax.legend(frameon=False, fontsize=9)

    # (b) fold-enrichment of genes and methylation along ori->ter
    ax = axes[0, 1]
    gc, gf = fold_enrichment(g_signed, nbins)
    mc, mf = fold_enrichment(m_signed, nbins)
    ax.plot(gc, gf, color="#2c3e50", lw=1.6, label="gene density")
    ax.plot(mc, mf, color="#c0392b", lw=1.6, label="methylation density")
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.axvline(0, color="#27ae60", lw=1.4)
    ax.text(0, ax.get_ylim()[1], "oriC", color="#27ae60", ha="center",
            va="bottom", fontsize=9)
    for xt in (-1, 1):
        ax.axvline(xt, color="#c0392b", lw=1.2, ls=":")
    ax.text(1, ax.get_ylim()[1], "ter", color="#c0392b", ha="right",
            va="bottom", fontsize=9)
    ax.set_xlabel("ori$\\rightarrow$ter coordinate (signed; replichores)")
    ax.set_ylabel("density / replichore mean")
    ax.set_title("(b) Spatial density vs replication axis")
    ax.legend(frameon=False, fontsize=9)

    # (c) methylation level vs distance from oriC
    ax = axes[1, 0]
    if m_absfrac:
        a = np.array([x for x, _ in m_absfrac])
        fr = np.array([f for _, f in m_absfrac])
        edges = np.linspace(0, 1, nbins // 2 + 1)
        idx = np.clip(np.digitize(a, edges) - 1, 0, len(edges) - 2)
        cen = (edges[:-1] + edges[1:]) / 2
        means = np.array([fr[idx == b].mean() if np.any(idx == b) else np.nan
                          for b in range(len(cen))])
        sems = np.array([fr[idx == b].std() / max(1, np.sqrt((idx == b).sum()))
                         if np.any(idx == b) else np.nan for b in range(len(cen))])
        ax.errorbar(cen, means, yerr=sems, fmt="o-", color="#8e44ad",
                    ms=4, lw=1.3, capsize=2)
        ax.set_ylim(0, 1)
    ax.set_xlabel("distance from oriC (0=oriC, 1=ter; replichores folded)")
    ax.set_ylabel("mean fraction modified")
    ax.set_title("(c) Methylation level gradient")

    # (d) gene co-orientation vs distance from oriC
    ax = axes[1, 1]
    if co_abs:
        a = np.array([x for x, _ in co_abs])
        co = np.array([1 if c else 0 for _, c in co_abs])
        edges = np.linspace(0, 1, nbins // 2 + 1)
        idx = np.clip(np.digitize(a, edges) - 1, 0, len(edges) - 2)
        cen = (edges[:-1] + edges[1:]) / 2
        frac = np.array([co[idx == b].mean() if np.any(idx == b) else np.nan
                         for b in range(len(cen))])
        ax.bar(cen, frac, width=(1.0 / (len(cen))) * 0.9, color="#16a085")
        ax.axhline(0.5, color="#999", ls="--", lw=0.8)
        ax.set_ylim(0, 1)
    ax.set_xlabel("distance from oriC (0=oriC, 1=ter; replichores folded)")
    ax.set_ylabel("fraction co-oriented (leading strand)")
    ax.set_title("(d) Gene orientation vs replication")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def make_circular_genome(outpath, genome, rep_order, headers, feats_by, calls_by,
                         orit_by, motif_occ, motif_list, binbp, plot_fmt,
                         include_motifs):
    """Circos-style whole-genome view. Concentric tracks (outer->inner):
    replicon ideogram + oriC/ter markers, gene density (+ strand), gene density
    (- strand), methylation density, GC skew (diverging). Optional motif-density
    rings are drawn outside the ideogram, one per motif."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import numpy as np

    rep_len = {k: len(genome[k]) for k in rep_order}
    total = sum(rep_len.values()) or 1
    rep_color = {k: plt.cm.tab10(i % 10) for i, k in enumerate(rep_order)}
    lab = lambda k: headers.get(k, k)  # noqa: E731

    gap = 2 * math.pi * 0.02
    usable = 2 * math.pi - len(rep_order) * gap
    layout = {}
    cur = math.pi / 2
    for k in rep_order:
        span = usable * rep_len[k] / total
        layout[k] = (cur, cur - span)
        cur -= span + gap

    def bars(positions, k):
        L = rep_len[k]
        t0, t1 = layout[k]
        nb = max(2, int(math.ceil(L / binbp)))
        counts, edges = np.histogram(positions, bins=nb, range=(0, L))
        te = t0 + (edges / L) * (t1 - t0)
        centers = (te[:-1] + te[1:]) / 2
        widths = np.abs(np.diff(te))
        return centers, widths, counts

    def skew_bins(k):
        seq = genome[k]
        L = rep_len[k]
        t0, t1 = layout[k]
        nb = max(2, int(math.ceil(L / binbp)))
        edges = np.linspace(0, L, nb + 1).astype(int)
        vals = []
        for i in range(nb):
            w = seq[edges[i]:edges[i + 1]]
            gc = w.count("G") + w.count("C")
            vals.append((w.count("G") - w.count("C")) / (gc if gc else 1))
        te = t0 + (edges.astype(float) / L) * (t1 - t0)
        centers = (te[:-1] + te[1:]) / 2
        widths = np.abs(np.diff(te))
        return centers, widths, np.array(vals)

    per = {}
    gene_max = meth_max = 1
    for k in rep_order:
        feats = feats_by.get(k, [])
        plus = [(f.start + f.end) // 2 for f in feats if f.strand == "+"]
        minus = [(f.start + f.end) // 2 for f in feats if f.strand == "-"]
        mpos = [c.pos for c in calls_by.get(k, [])]
        per[k] = dict(plus=bars(plus, k), minus=bars(minus, k),
                      meth=bars(mpos, k), skew=skew_bins(k))
        gene_max = max(gene_max, per[k]["plus"][2].max() if per[k]["plus"][2].size else 0,
                       per[k]["minus"][2].max() if per[k]["minus"][2].size else 0)
        meth_max = max(meth_max, per[k]["meth"][2].max() if per[k]["meth"][2].size else 0)
    skew_absmax = max((np.abs(per[k]["skew"][2]).max() for k in rep_order), default=1) or 1

    motif_used = motif_list if (include_motifs and motif_list) else []
    motif_bars, motif_max = {m: {} for m in motif_used}, {}
    for m in motif_used:
        mm = 1
        for k in rep_order:
            cb = bars(motif_occ.get(m, {}).get(k, []), k)
            motif_bars[m][k] = cb
            mm = max(mm, cb[2].max() if cb[2].size else 0)
        motif_max[m] = mm

    fig = plt.figure(figsize=(10.5, 10.5))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.axis("off")

    R_GC_BASE, GC_AMP = 0.32, 0.12
    R_METH, H_METH = 0.48, 0.18
    R_GMINUS, H_G = 0.70, 0.12
    R_GPLUS = 0.84
    R_IDEO, W_IDEO = 0.98, 0.05
    R_MARK = R_IDEO + W_IDEO + 0.01
    R_MOTIF0 = R_IDEO + W_IDEO + 0.09
    motif_h = 0.085
    outer = R_MOTIF0 + len(motif_used) * (motif_h + 0.02) + 0.10

    theta = np.linspace(0, 2 * math.pi, 400)
    ax.plot(theta, [R_GC_BASE] * len(theta), color="#cccccc", lw=0.5, zorder=1)

    for k in rep_order:
        t0, t1 = layout[k]
        tc = (t0 + t1) / 2
        ax.bar(tc, W_IDEO, width=abs(t1 - t0), bottom=R_IDEO, color=rep_color[k],
               edgecolor="white", linewidth=0.6, align="center", zorder=5)
        # GC skew (diverging from baseline)
        c, w, vals = per[k]["skew"]
        ax.bar(c, (np.clip(vals, 0, None) / skew_absmax) * GC_AMP, width=w,
               bottom=R_GC_BASE, color="#2ecc71", linewidth=0, align="center", zorder=2)
        negh = (np.clip(-vals, 0, None) / skew_absmax) * GC_AMP
        ax.bar(c, negh, width=w, bottom=R_GC_BASE - negh, color="#8e44ad",
               linewidth=0, align="center", zorder=2)
        # methylation density
        c, w, cnt = per[k]["meth"]
        ax.bar(c, (cnt / meth_max) * H_METH, width=w, bottom=R_METH,
               color="#c0392b", linewidth=0, align="center", zorder=2)
        # genes (- strand) then (+ strand)
        c, w, cnt = per[k]["minus"]
        ax.bar(c, (cnt / gene_max) * H_G, width=w, bottom=R_GMINUS,
               color="#95a5a6", linewidth=0, align="center", zorder=2)
        c, w, cnt = per[k]["plus"]
        ax.bar(c, (cnt / gene_max) * H_G, width=w, bottom=R_GPLUS,
               color="#2c3e50", linewidth=0, align="center", zorder=2)
        # oriC/ter markers
        oric, ter = orit_by[k]
        for posm, colm, labm in ((oric, "#27ae60", "ori"), (ter, "#c0392b", "ter")):
            th = t0 + (posm / rep_len[k]) * (t1 - t0)
            ax.plot([th, th], [R_IDEO - 0.015, R_MARK + 0.03], color=colm, lw=2.2,
                    zorder=6)
            ax.text(th, R_MARK + 0.06, labm, color=colm, ha="center", va="center",
                    fontsize=7, fontweight="bold", zorder=6)
        # motif rings (outside ideogram)
        for mi, m in enumerate(motif_used):
            rb = R_MOTIF0 + mi * (motif_h + 0.02)
            c, w, cnt = motif_bars[m][k]
            ax.bar(c, (cnt / motif_max[m]) * motif_h, width=w, bottom=rb,
                   color=plt.cm.Set2(mi % 8), linewidth=0, align="center", zorder=2)
        # replicon label
        deg = math.degrees(tc) % 360
        ha = "left" if (deg < 90 or deg > 270) else "right"
        ax.text(tc, outer + 0.02, f"{lab(k)}\n{rep_len[k] / 1e6:.2f} Mb", ha=ha,
                va="center", fontsize=9, color=rep_color[k], fontweight="bold")

    ax.set_ylim(0, outer + 0.12)
    legend = [Patch(facecolor="#2c3e50", label="genes (+ strand)"),
              Patch(facecolor="#95a5a6", label="genes ($-$ strand)"),
              Patch(facecolor="#c0392b", label="methylation density"),
              Patch(facecolor="#2ecc71", label="GC skew (+)"),
              Patch(facecolor="#8e44ad", label="GC skew ($-$)")]
    for mi, m in enumerate(motif_used):
        legend.append(Patch(facecolor=plt.cm.Set2(mi % 8), label=f"motif {m}"))
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=3, frameon=False, fontsize=8)
    ax.set_title("Genome-wide distribution of genes, methylation"
                 + (", motifs" if motif_used else "") + "  (oriC/ter marked)",
                 fontsize=13, fontweight="bold", pad=24)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def make_motif_figure(outpath, repname, motif_signed, nbins, plot_fmt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.6))
    cmap = plt.cm.tab10
    for i, (m, vals) in enumerate(motif_signed.items()):
        if not vals:
            continue
        c, f = fold_enrichment(vals, nbins)
        ax.plot(c, f, lw=1.6, color=cmap(i % 10), label=f"{m} (n={len(vals)})")
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.axvline(0, color="#27ae60", lw=1.2)
    ax.set_xlabel("ori$\\rightarrow$ter coordinate (signed)")
    ax.set_ylabel("methylated-on-motif density / mean")
    ax.set_title(f"{repname}: methylation-on-motif spatial distribution")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def parse_locus_overrides(s):
    out = {}
    if not s:
        return out
    for tok in s.split(","):
        tok = tok.strip()
        if ":" in tok:
            k, v = tok.rsplit(":", 1)
            out[k.strip()] = int(v)
    return out


def replicon_label(genome_headers, seqid):
    return genome_headers.get(seqid, seqid)


def read_headers(path, strip):
    labels = {}
    if not path:
        return labels
    with smart_open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            h = line[1:].strip()
            sid = h.split()[0]
            if strip:
                sid = strip_version(sid)
            d = h.lower()
            lab = sid
            if "chromosome" in d:
                lab = "chromosome"
            elif "plasmid" in d:
                after = d.split("plasmid", 1)[1].strip()
                tok = after.split(",")[0].split()[0] if after else ""
                lab = tok.upper() if tok else "plasmid"
            labels[sid] = lab
    return labels


def build_argparser():
    p = argparse.ArgumentParser(
        description="Spatial distribution of methylation and features vs oriC/ter.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--gff", required=True)
    p.add_argument("--fasta", required=True,
                   help="Reference FASTA (required: oriC/ter come from GC skew).")
    p.add_argument("--methylation", required=True)
    p.add_argument("--meth-format", default="auto",
                   choices=["auto", "bedmethyl", "gff"])
    p.add_argument("--motifs", default=None,
                   help="Optional motif list; adds a per-motif methylation track.")
    p.add_argument("--mod-codes", default=None,
                   help="Comma-separated modification codes to keep (e.g. 'a').")
    p.add_argument("--feature-types", default=None,
                   help="Comma-separated GFF types to keep (e.g. 'CDS').")
    p.add_argument("--id-attr", default="ID")
    p.add_argument("--min-coverage", type=float, default=0.0)
    p.add_argument("--min-frac", type=float, default=0.5)
    p.add_argument("--skew-window", type=int, default=10000,
                   help="Window (bp) for cumulative GC-skew oriC/ter detection.")
    p.add_argument("--oric", default=None,
                   help="Override oriC, format 'seqid:pos[,seqid:pos]' (0-based).")
    p.add_argument("--ter", default=None,
                   help="Override ter, format 'seqid:pos[,seqid:pos]' (0-based).")
    p.add_argument("--bins", type=int, default=40,
                   help="Bins across the signed ori->ter axis [-1,1].")
    p.add_argument("--min-skew-amp", type=int, default=500,
                   help="Warn if cumulative skew amplitude is below this (weak "
                        "oriC/ter signal).")
    p.add_argument("--strip-version", action="store_true",
                   help="Strip trailing .N from seqids (auto-applied if it "
                        "reconciles GFF and FASTA).")
    p.add_argument("--plots", action="store_true")
    p.add_argument("--circular", action="store_true",
                   help="Also render the whole-genome circular plot (implied by "
                        "--plots).")
    p.add_argument("--no-motif-rings", action="store_true",
                   help="Do not add motif-density rings to the circular plot even "
                        "if --motifs is given.")
    p.add_argument("--plot-format", default="png", choices=["png", "pdf", "svg"])
    p.add_argument("--outdir", default="oriter_out")
    p.add_argument("--prefix", default="sample")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv=None):
    args = build_argparser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    st = try_scipy()
    feature_types = set(args.feature_types.split(",")) if args.feature_types else None
    keep_codes = set(args.mod_codes.split(",")) if args.mod_codes else None

    sys.stderr.write("Parsing GFF / FASTA / methylation ...\n")
    features, embedded = parse_gff3(args.gff, feature_types, args.id_attr)
    genome = load_fasta(args.fasta) if args.fasta else embedded
    calls = parse_methylation(args.methylation, args.meth_format, keep_codes)

    # seqid reconciliation (auto strip a pure version-suffix mismatch)
    gff_ids = {f.seqid for f in features}
    fa_ids = set(genome)
    apply_strip = args.strip_version
    if not apply_strip and gff_ids and not (gff_ids & fa_ids):
        if {strip_version(x) for x in gff_ids} & {strip_version(x) for x in fa_ids}:
            sys.stderr.write("[auto] reconciling seqids by stripping '.N' suffix.\n")
            apply_strip = True
    if apply_strip:
        for f in features:
            f.seqid = strip_version(f.seqid)
        for c in calls:
            c.seqid = strip_version(c.seqid)
        genome = {strip_version(k): v for k, v in genome.items()}

    headers = read_headers(args.fasta, apply_strip)
    if not (set(genome) & {f.seqid for f in features}):
        sys.stderr.write("[ERROR] GFF and FASTA seqids do not match; cannot map "
                         "features. Try --strip-version.\n")
        return 3

    oric_ov = parse_locus_overrides(args.oric)
    ter_ov = parse_locus_overrides(args.ter)

    sys.stderr.write(f"  {len(features)} features, {len(calls)} calls, "
                     f"{len(genome)} replicon(s).\n")

    # optional motif occurrences (positions for rings) + covered sets (membership)
    motif_list = parse_motifs(args.motifs) if args.motifs else []
    motif_occ, motif_cov = (find_motif_occurrences(genome, motif_list)
                            if motif_list else ({}, {}))

    # group features / calls by replicon
    feats_by = defaultdict(list)
    for f in features:
        feats_by[f.seqid].append(f)
    calls_by = defaultdict(list)
    for c in calls:
        if (c.coverage is None) or (c.coverage >= args.min_coverage):
            calls_by[c.seqid].append(c)

    rep_order = sorted(genome, key=lambda x: -len(genome[x]))

    stats_rows = []
    bin_rows = []
    # genome-wide accumulators for co-orientation + KS
    G_co_k = G_co_n = 0
    G_meth_abs = []
    G_gene_abs = []
    G_frac_pairs = []
    orit_by = {}  # seqid -> (oriC, ter) for the circular plot

    for seqid in rep_order:
        seq = genome[seqid]
        L = len(seq)
        label = replicon_label(headers, seqid)
        sp = gc_skew(seq, args.skew_window)
        oric = oric_ov.get(seqid, oric_ov.get(strip_version(seqid), sp.oric))
        ter = ter_ov.get(seqid, ter_ov.get(strip_version(seqid), sp.ter))
        # rebuild SkewProfile with possibly-overridden markers for plotting
        sp = SkewProfile(sp.centers, sp.cum, oric, ter, sp.amplitude)
        orit_by[seqid] = (oric, ter)
        weak = sp.amplitude < args.min_skew_amp
        if weak:
            sys.stderr.write(f"[warn] {label}: weak GC-skew (amp={sp.amplitude}); "
                             f"oriC/ter may be unreliable.\n")

        # ---- genes ----
        g_signed, g_abs = [], []
        co_abs = []
        co_k = co_n = 0
        for f in feats_by.get(seqid, []):
            mid = (f.start + f.end) // 2
            s, a, arm = ori_ter(mid, oric, ter, L)
            g_signed.append(s)
            g_abs.append(a)
            cob = co_oriented(f.strand, arm)
            if cob is not None:
                co_n += 1
                co_k += 1 if cob else 0
                co_abs.append((a, cob))
        G_co_k += co_k
        G_co_n += co_n
        G_gene_abs += g_abs

        # ---- methylation ----
        m_signed, m_abs = [], []
        m_absfrac = []
        meth_considered = calls_by.get(seqid, [])
        for c in meth_considered:
            s, a, arm = ori_ter(c.pos, oric, ter, L)
            m_signed.append(s)
            m_abs.append(a)
            if c.frac is not None:
                m_absfrac.append((a, c.frac))
        G_meth_abs += m_abs
        G_frac_pairs += m_absfrac

        # ---- statistics ----
        co_frac = (co_k / co_n) if co_n else float("nan")
        co_p = float("nan")
        if st and co_n:
            co_p = st.binomtest(co_k, co_n, 0.5).pvalue
        gene_ks_p = meth_ks_p = float("nan")
        gene_mean = (sum(g_abs) / len(g_abs)) if g_abs else float("nan")
        meth_mean = (sum(m_abs) / len(m_abs)) if m_abs else float("nan")
        if st:
            if len(g_abs) > 2:
                gene_ks_p = st.kstest(g_abs, "uniform").pvalue
            if len(m_abs) > 2:
                meth_ks_p = st.kstest(m_abs, "uniform").pvalue
        level_rho = level_p = float("nan")
        if st and len(m_absfrac) > 10:
            xa = [x for x, _ in m_absfrac]
            fa = [f for _, f in m_absfrac]
            level_rho, level_p = st.spearmanr(xa, fa)

        stats_rows.append(dict(
            replicon=seqid, label=label, length=L,
            oriC=oric, ter=ter, skew_amplitude=sp.amplitude, weak_skew=int(weak),
            n_features=len(feats_by.get(seqid, [])),
            n_meth_calls=len(meth_considered),
            gene_co_oriented_frac=co_frac, gene_co_oriented_n=co_n,
            gene_co_binom_p=co_p,
            gene_mean_oridist=gene_mean, gene_uniform_KS_p=gene_ks_p,
            meth_mean_oridist=meth_mean, meth_uniform_KS_p=meth_ks_p,
            meth_level_spearman_rho=level_rho, meth_level_spearman_p=level_p,
        ))

        # ---- binned profile dump ----
        try:
            import numpy as np
            nb = args.bins
            edges = np.linspace(-1, 1, nb + 1)
            gc_counts, _ = np.histogram(g_signed, bins=edges)
            mc_counts, _ = np.histogram(m_signed, bins=edges)
            cen = (edges[:-1] + edges[1:]) / 2
            # mean frac and co-orientation per |bin| done on folded axis separately
            for b in range(nb):
                bin_rows.append(dict(replicon=seqid, label=label,
                                     signed_center=round(float(cen[b]), 4),
                                     n_genes=int(gc_counts[b]),
                                     n_meth=int(mc_counts[b])))
        except ImportError:
            pass

        # ---- plots ----
        if args.plots:
            outpath = os.path.join(
                args.outdir, f"{args.prefix}.{label}.spatial.{args.plot_format}")
            make_replicon_figure(outpath, label, L, sp, g_signed, m_signed,
                                 m_absfrac, co_abs, args.bins, args.plot_format)
            sys.stderr.write(f"  wrote {outpath}\n")

            if motif_list:
                motif_signed = {m: [] for m in motif_list}
                for c in meth_considered:
                    if c.frac is not None and c.frac < args.min_frac:
                        continue
                    for m in motif_list:
                        if (c.seqid, c.pos) in motif_cov[m]:
                            s, _a, _arm = ori_ter(c.pos, oric, ter, L)
                            motif_signed[m].append(s)
                mpath = os.path.join(
                    args.outdir,
                    f"{args.prefix}.{label}.motif_spatial.{args.plot_format}")
                make_motif_figure(mpath, label, motif_signed, args.bins,
                                  args.plot_format)
                sys.stderr.write(f"  wrote {mpath}\n")

    # ---- whole-genome circular plot ----
    if args.plots or args.circular:
        binbp = max(2000, sum(len(genome[k]) for k in genome) // 1500)
        cpath = os.path.join(args.outdir,
                             f"{args.prefix}.circular_genome.{args.plot_format}")
        make_circular_genome(cpath, genome, rep_order, headers, feats_by, calls_by,
                             orit_by, motif_occ, motif_list, binbp, args.plot_format,
                             include_motifs=(bool(motif_list)
                                             and not args.no_motif_rings))
        sys.stderr.write(f"  wrote {cpath}\n")

    # genome-wide row
    g_co_frac = (G_co_k / G_co_n) if G_co_n else float("nan")
    g_co_p = st.binomtest(G_co_k, G_co_n, 0.5).pvalue if (st and G_co_n) else float("nan")
    g_gene_ks = st.kstest(G_gene_abs, "uniform").pvalue if (st and len(G_gene_abs) > 2) else float("nan")
    g_meth_ks = st.kstest(G_meth_abs, "uniform").pvalue if (st and len(G_meth_abs) > 2) else float("nan")
    g_rho = g_rp = float("nan")
    if st and len(G_frac_pairs) > 10:
        g_rho, g_rp = st.spearmanr([x for x, _ in G_frac_pairs],
                                   [f for _, f in G_frac_pairs])
    stats_rows.append(dict(
        replicon="GENOME", label="all replicons",
        length=sum(len(genome[k]) for k in genome),
        oriC="", ter="", skew_amplitude="", weak_skew="",
        n_features=len(features), n_meth_calls=sum(len(v) for v in calls_by.values()),
        gene_co_oriented_frac=g_co_frac, gene_co_oriented_n=G_co_n,
        gene_co_binom_p=g_co_p,
        gene_mean_oridist=(sum(G_gene_abs) / len(G_gene_abs)) if G_gene_abs else float("nan"),
        gene_uniform_KS_p=g_gene_ks,
        meth_mean_oridist=(sum(G_meth_abs) / len(G_meth_abs)) if G_meth_abs else float("nan"),
        meth_uniform_KS_p=g_meth_ks,
        meth_level_spearman_rho=g_rho, meth_level_spearman_p=g_rp))

    # write stats
    spath = os.path.join(args.outdir, f"{args.prefix}.oriter_stats.tsv")
    cols = ["replicon", "label", "length", "oriC", "ter", "skew_amplitude",
            "weak_skew", "n_features", "n_meth_calls",
            "gene_co_oriented_frac", "gene_co_oriented_n", "gene_co_binom_p",
            "gene_mean_oridist", "gene_uniform_KS_p",
            "meth_mean_oridist", "meth_uniform_KS_p",
            "meth_level_spearman_rho", "meth_level_spearman_p"]
    with open(spath, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        wr.writeheader()
        for r in stats_rows:
            wr.writerow({k: (f"{r[k]:.4g}" if isinstance(r[k], float) else r[k])
                         for k in cols})
    bpath = os.path.join(args.outdir, f"{args.prefix}.spatial_bins.tsv")
    with open(bpath, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["replicon", "label", "signed_center",
                                            "n_genes", "n_meth"], delimiter="\t")
        wr.writeheader()
        for r in bin_rows:
            wr.writerow(r)

    # console summary
    sys.stderr.write("\n=== oriC/ter spatial summary ===\n")
    for r in stats_rows:
        if r["replicon"] == "GENOME":
            continue
        sys.stderr.write(
            f"  {r['label']:<12} oriC={r['oriC']:>9} ter={r['ter']:>9} | "
            f"co-oriented {r['gene_co_oriented_frac']*100:.1f}% "
            f"(p={r['gene_co_binom_p']:.1e}) | "
            f"meth level vs ori rho={r['meth_level_spearman_rho']:+.3f} "
            f"(p={r['meth_level_spearman_p']:.1e})\n")
    sys.stderr.write(f"  stats -> {spath}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
