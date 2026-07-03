#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
motif_methylation_cross.py

Cross-reference DNA methylation basecalls, user-supplied DNA motifs, and a
genomic annotation (GFF3) to ask, per annotated feature, whether observed
methylation is associated with the presence of a given motif.

WHAT IT DOES
------------
1. Parses a GFF3 annotation (feature coordinates, strand, type, attributes).
2. Loads the reference sequence (from --fasta, or from a ``##FASTA`` block
   embedded in the GFF3) -- required to locate motifs.
3. Reads a newline-delimited list of DNA motifs (IUPAC allowed). Each motif may
   optionally carry the 0-based offset of the modified base within the motif
   (e.g. ``GATC<TAB>1`` for Dam: the A at index 1). Lines beginning with '#'
   are ignored.
4. Performs naive (overlap-aware) matching of each motif on BOTH strands of
   every replicon, recording for each occurrence: replicon, coordinates,
   strand, and -- when an offset is given -- the exact genomic position of the
   modified base.
5. Parses the methylation file:
       * bedMethyl (modkit-style, default for .bed)   -> has strand, coverage, %mod
       * GFF3-style modification calls (default for .gff/.gff3)
6. Crosses the three layers, taking REPLICON, STRAND and FEATURE into account,
   and reports for each feature:
       * methylation calls falling inside the feature,
       * how many of those fall on a motif (and, if offsets are given, on the
         exact modified base of the motif),
       * whether a methylated feature carries its methylation on the motif.

WHAT IT IS (and is NOT)
-----------------------
The matching is intentionally *naive*: a methylation call is "on a motif" if its
genomic position (strand-aware by default) lies within a motif occurrence. This
is a descriptive cross-tabulation. The optional Fisher tests quantify
*association*, not causation -- assay/coverage bias, motif palindromy and
annotation overlap all shape the result and are reported so they can be judged.

OUTPUTS (written to --outdir, prefixed with --prefix)
-----------------------------------------------------
  <prefix>.methylation_calls.tsv  per-call annotation
  <prefix>.per_feature.tsv        per-feature aggregate
  <prefix>.per_feature_motif.tsv  per-(feature, motif) aggregate
  <prefix>.motif_summary.tsv      per-motif genome-wide aggregate
  <prefix>.replicons.tsv          replicon lengths (for re-plotting)
  <prefix>.motif_hits.tsv         every motif occurrence (replicon/strand/coords)
  <prefix>.enrichment.tsv         Fisher tests (only with --fisher)
  With --plots (requires matplotlib):
  <prefix>.circular.<motif>.png            motif occurrences + methylated-on-motif
                                           density over the genome, per motif
  <prefix>.motif_methylation_summary.png   per-motif rate, site enrichment, on/off
  <prefix>.feature_motif_association.png   P(feature methylated | motif present)

Sequence-id reconciliation: motifs are located on the FASTA coordinates, so GFF
and methylation must share the FASTA seqids. A pure trailing '.N' version-suffix
mismatch (e.g. GFF 'X' vs FASTA 'X.1') is auto-resolved; --strip-version forces
it; any other mismatch is reported and aborts rather than silently yielding
all-zero per-feature counts.

Author: built for prokaryotic ONT methylation / annotation workflows.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

__version__ = "1.2.0"


# --------------------------------------------------------------------------- #
# IUPAC handling
# --------------------------------------------------------------------------- #
IUPAC_TO_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]",
}

IUPAC_COMPLEMENT = {
    "A": "T", "T": "A", "U": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
    "N": "N",
}


def reverse_complement_iupac(motif: str) -> str:
    """Reverse-complement an IUPAC motif (used to find minus-strand hits)."""
    try:
        return "".join(IUPAC_COMPLEMENT[b] for b in reversed(motif.upper()))
    except KeyError as exc:  # pragma: no cover - defensive
        raise ValueError(f"Illegal IUPAC base in motif {motif!r}: {exc}") from exc


def motif_to_regex(motif: str) -> "re.Pattern":
    """
    Compile an IUPAC motif into an overlap-aware regex.

    A lookahead is used so that overlapping occurrences are all reported
    (e.g. 'AA' in 'AAA' yields two hits). Motifs are fixed-length, so the match
    span equals len(motif).
    """
    body = "".join(IUPAC_TO_REGEX[b] for b in motif.upper())
    return re.compile(f"(?=({body}))")


# --------------------------------------------------------------------------- #
# Data models. All internal coordinates are 0-based, inclusive [start, end].
# --------------------------------------------------------------------------- #
@dataclass
class Feature:
    seqid: str
    ftype: str
    start: int          # 0-based inclusive
    end: int            # 0-based inclusive
    strand: str         # '+', '-' or '.'
    fid: str
    attrs: Dict[str, str] = field(default_factory=dict)


@dataclass
class MotifHit:
    seqid: str
    start: int          # 0-based inclusive
    end: int            # 0-based inclusive
    strand: str         # '+' or '-'
    motif: str          # motif as written (5'->3')
    mod_pos: Optional[int]  # 0-based genomic position of modified base, if known


@dataclass
class MethCall:
    seqid: str
    pos: int            # 0-based position of the modified base
    strand: str         # '+', '-' or '.'
    mod_code: str
    coverage: Optional[float]
    frac: Optional[float]   # fraction modified in [0, 1], or None if unknown


# --------------------------------------------------------------------------- #
# IO helpers
# --------------------------------------------------------------------------- #
def smart_open(path: str):
    """Open plain or gzipped text transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


# --------------------------------------------------------------------------- #
# GFF3 parsing (features + optional embedded FASTA)
# --------------------------------------------------------------------------- #
def parse_gff3(
    path: str,
    feature_types: Optional[set],
    id_attr: str,
) -> Tuple[List[Feature], Dict[str, str]]:
    """
    Parse a GFF3 file.

    Returns (features, embedded_fasta) where embedded_fasta maps seqid->sequence
    for any sequence found after a ``##FASTA`` directive (may be empty).
    """
    features: List[Feature] = []
    embedded_fasta: Dict[str, str] = {}
    in_fasta = False
    cur_id: Optional[str] = None
    cur_seq: List[str] = []
    auto_n = 0

    def parse_attrs(field9: str) -> Dict[str, str]:
        attrs: Dict[str, str] = {}
        for kv in field9.strip().split(";"):
            kv = kv.strip()
            if not kv:
                continue
            if "=" in kv:
                k, v = kv.split("=", 1)
            elif " " in kv:  # tolerate GTF-ish "key value"
                k, v = kv.split(" ", 1)
                v = v.strip('"')
            else:
                continue
            attrs[k.strip()] = v.strip()
        return attrs

    with smart_open(path) as fh:
        for line in fh:
            if in_fasta:
                if line.startswith(">"):
                    if cur_id is not None:
                        embedded_fasta[cur_id] = "".join(cur_seq).upper()
                    cur_id = line[1:].strip().split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line.strip())
                continue

            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attr_field = cols[:9]
            if feature_types and ftype not in feature_types:
                continue
            try:
                start0 = int(start) - 1
                end0 = int(end) - 1
            except ValueError:
                continue
            attrs = parse_attrs(attr_field)
            fid = attrs.get(id_attr)
            if fid is None:
                fid = attrs.get("locus_tag") or attrs.get("Name")
            if fid is None:
                auto_n += 1
                fid = f"{ftype}_{seqid}_{start0+1}_{end0+1}_{auto_n}"
            features.append(
                Feature(seqid=seqid, ftype=ftype, start=start0, end=end0,
                        strand=strand, fid=fid, attrs=attrs)
            )

    if in_fasta and cur_id is not None:
        embedded_fasta[cur_id] = "".join(cur_seq).upper()

    return features, embedded_fasta


# --------------------------------------------------------------------------- #
# FASTA loading
# --------------------------------------------------------------------------- #
def load_fasta(path: str) -> Dict[str, str]:
    genome: Dict[str, str] = {}
    cur_id: Optional[str] = None
    cur_seq: List[str] = []
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    genome[cur_id] = "".join(cur_seq).upper()
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_id is not None:
        genome[cur_id] = "".join(cur_seq).upper()
    return genome


# --------------------------------------------------------------------------- #
# Motif file parsing
# --------------------------------------------------------------------------- #
def parse_motifs(path: str) -> List[Tuple[str, Optional[int]]]:
    """
    One motif per line. Optional second whitespace-separated field = 0-based
    offset of the modified base within the motif. Lines starting with '#'
    are ignored.
    """
    motifs: List[Tuple[str, Optional[int]]] = []
    with smart_open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            motif = parts[0].upper()
            for b in motif:
                if b not in IUPAC_TO_REGEX:
                    raise ValueError(f"Illegal IUPAC base {b!r} in motif {motif!r}")
            offset: Optional[int] = None
            if len(parts) > 1:
                try:
                    offset = int(parts[1])
                    if not (0 <= offset < len(motif)):
                        raise ValueError
                except ValueError:
                    raise ValueError(
                        f"Bad modified-base offset for motif {motif!r}: {parts[1]!r} "
                        f"(must be an integer in [0, {len(motif)-1}])"
                    )
            motifs.append((motif, offset))
    if not motifs:
        raise ValueError("No motifs parsed from motif file.")
    return motifs


# --------------------------------------------------------------------------- #
# Motif search on both strands
# --------------------------------------------------------------------------- #
def find_motif_hits(
    genome: Dict[str, str],
    motifs: List[Tuple[str, Optional[int]]],
) -> List[MotifHit]:
    """
    Naive overlap-aware search of every motif on both strands of every replicon.

    Plus strand: search the motif directly; modified base at start + offset.
    Minus strand: search the reverse complement of the motif in the plus
    sequence. An occurrence at plus-coords [m_start, m_end] is a minus-strand
    motif; because revcomp reverses orientation, the motif's 5' base maps to
    m_end, so the modified base maps to m_end - offset.
    """
    hits: List[MotifHit] = []
    for motif, offset in motifs:
        L = len(motif)
        fwd_re = motif_to_regex(motif)
        rc = reverse_complement_iupac(motif)
        is_palindrome = (rc == motif)
        rev_re = fwd_re if is_palindrome else motif_to_regex(rc)
        for seqid, seq in genome.items():
            # plus strand
            for m in fwd_re.finditer(seq):
                s = m.start()
                mod = (s + offset) if offset is not None else None
                hits.append(MotifHit(seqid, s, s + L - 1, "+", motif, mod))
            # minus strand
            for m in rev_re.finditer(seq):
                s = m.start()
                e = s + L - 1
                mod = (e - offset) if offset is not None else None
                hits.append(MotifHit(seqid, s, e, "-", motif, mod))
    return hits


# --------------------------------------------------------------------------- #
# Methylation parsing
# --------------------------------------------------------------------------- #
def _norm_frac(value: float) -> float:
    """Normalise a fraction that may be given as a percentage."""
    return value / 100.0 if value > 1.0 else value


def parse_bedmethyl(
    path: str,
    keep_codes: Optional[set],
) -> List[MethCall]:
    """
    Parse a modkit-style bedMethyl file.

    Column layout (0-based index):
      0 chrom | 1 start | 2 end | 3 mod_code | 4 score | 5 strand
      9 Nvalid_cov | 10 fraction_modified(%) ...
    A 6-column generic BED is tolerated (no coverage / fraction -> presence).
    """
    calls: List[MethCall] = []
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                c = line.split()
            if len(c) < 6:
                continue
            chrom = c[0]
            try:
                start = int(c[1])
            except ValueError:
                continue
            strand = c[5] if c[5] in ("+", "-", ".") else "."
            mod_code = c[3] if len(c) > 3 else "."
            if keep_codes and mod_code not in keep_codes:
                continue
            coverage: Optional[float] = None
            frac: Optional[float] = None
            if len(c) >= 11:
                try:
                    coverage = float(c[9])
                except ValueError:
                    coverage = None
                try:
                    frac = _norm_frac(float(c[10]))
                except ValueError:
                    frac = None
            calls.append(MethCall(chrom, start, strand, mod_code, coverage, frac))
    return calls


def parse_gff_meth(
    path: str,
    keep_codes: Optional[set],
) -> List[MethCall]:
    """
    Parse GFF3-style modification calls. The feature 'type' is treated as the
    modification code. Coverage / fraction are read from attributes when present
    (keys tried: coverage/cov/Nvalid_cov ; frac/frac_modified/fraction/
    percent_modified). IPDRatio etc. are left untouched.
    """
    cov_keys = ("coverage", "cov", "Nvalid_cov", "valid_coverage", "Ncov")
    frac_keys = ("frac", "frac_modified", "fraction", "percent_modified",
                 "fraction_modified")
    calls: List[MethCall] = []
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            chrom, _src, mtype, start, _end, _score, strand = c[:7]
            if keep_codes and mtype not in keep_codes:
                continue
            try:
                pos = int(start) - 1  # GFF is 1-based
            except ValueError:
                continue
            strand = strand if strand in ("+", "-", ".") else "."
            attrs = {}
            if len(c) >= 9:
                for kv in c[8].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        attrs[k.strip()] = v.strip()
            coverage = None
            for k in cov_keys:
                if k in attrs:
                    try:
                        coverage = float(attrs[k]); break
                    except ValueError:
                        pass
            frac = None
            for k in frac_keys:
                if k in attrs:
                    try:
                        frac = _norm_frac(float(attrs[k])); break
                    except ValueError:
                        pass
            calls.append(MethCall(chrom, pos, strand, mtype, coverage, frac))
    return calls


def parse_methylation(path: str, fmt: str, keep_codes: Optional[set]) -> List[MethCall]:
    if fmt == "auto":
        low = path.lower()
        if low.endswith((".bed", ".bed.gz", ".bedmethyl", ".bedmethyl.gz")):
            fmt = "bedmethyl"
        elif low.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz", ".gff2")):
            fmt = "gff"
        else:
            raise ValueError(
                f"Cannot infer methylation format from {path!r}; "
                f"pass --meth-format bedmethyl|gff."
            )
    if fmt == "bedmethyl":
        return parse_bedmethyl(path, keep_codes)
    if fmt == "gff":
        return parse_gff_meth(path, keep_codes)
    raise ValueError(f"Unknown methylation format: {fmt}")


# --------------------------------------------------------------------------- #
# Interval overlap by sweep line (point in interval, with multiplicity)
# --------------------------------------------------------------------------- #
def sweep_points_in_intervals(
    points: List[Tuple[int, int]],          # (pos, point_global_idx)
    intervals: List[Tuple[int, int, int]],  # (start, end, interval_global_idx)
) -> Dict[int, List[int]]:
    """
    For each point, return the list of interval indices that contain it.
    Coordinates are 0-based inclusive. O((P+I) log I + total_overlaps).
    """
    result: Dict[int, List[int]] = {}
    pts = sorted(points)
    ivs = sorted(intervals)
    active: List[Tuple[int, int]] = []  # heap of (end, idx)
    i, n = 0, len(ivs)
    for pos, pidx in pts:
        while i < n and ivs[i][0] <= pos:
            s, e, idx = ivs[i]
            heapq.heappush(active, (e, idx))
            i += 1
        while active and active[0][0] < pos:
            heapq.heappop(active)
        # remaining active intervals all satisfy start<=pos and end>=pos
        result[pidx] = [idx for (_e, idx) in active]
    return result


# --------------------------------------------------------------------------- #
# BH / FDR
# --------------------------------------------------------------------------- #
def benjamini_hochberg(pvals: List[float]) -> List[float]:
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        idx = order[rank]
        val = pvals[idx] * n / (rank + 1)
        prev = min(prev, val)
        q[idx] = min(prev, 1.0)
    return q


# --------------------------------------------------------------------------- #
# Core analysis
# --------------------------------------------------------------------------- #
@dataclass
class CallAnnotation:
    feature_idxs: List[int]
    hit_idxs: List[int]            # indices into the motif-hit list this call sits on
    motifs_hit: List[str]          # motif strings whose occurrence contains the call
    on_motif: bool                 # within any motif occurrence (strand-aware)
    on_mod_base: bool              # exactly on a motif modified base (needs offsets)
    considered: bool               # passes coverage filter
    methylated: bool               # passes fraction threshold (or presence)


def annotate_calls(
    calls: List[MethCall],
    features: List[Feature],
    hits: List[MotifHit],
    min_coverage: float,
    min_frac: float,
    ignore_motif_strand: bool,
    offsets_known: bool,
) -> List[CallAnnotation]:
    # ---- call -> features (strand-agnostic: genomic span containment) ----
    pts_by_seq: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for ci, call in enumerate(calls):
        pts_by_seq[call.seqid].append((call.pos, ci))
    feat_by_seq: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for fi, f in enumerate(features):
        feat_by_seq[f.seqid].append((f.start, f.end, fi))
    call_features: Dict[int, List[int]] = {}
    for seqid, pts in pts_by_seq.items():
        call_features.update(sweep_points_in_intervals(pts, feat_by_seq.get(seqid, [])))

    # ---- call -> motif occurrences (strand-aware unless ignore) ----
    def keyer(seqid: str, strand: str) -> Tuple[str, str]:
        return (seqid, "*") if ignore_motif_strand else (seqid, strand)

    pts_by_key: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    for ci, call in enumerate(calls):
        pts_by_key[keyer(call.seqid, call.strand)].append((call.pos, ci))
    hit_by_key: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = defaultdict(list)
    for hi, h in enumerate(hits):
        hit_by_key[keyer(h.seqid, h.strand)].append((h.start, h.end, hi))
    call_hits: Dict[int, List[int]] = {}
    for key, pts in pts_by_key.items():
        call_hits.update(sweep_points_in_intervals(pts, hit_by_key.get(key, [])))

    # ---- exact modified-base lookup (strand-aware) ----
    modbase: Dict[Tuple[str, str, int], set] = defaultdict(set)
    if offsets_known:
        for h in hits:
            if h.mod_pos is not None:
                k = (h.seqid, "*" if ignore_motif_strand else h.strand, h.mod_pos)
                modbase[k].add(h.motif)

    annotations: List[CallAnnotation] = []
    for ci, call in enumerate(calls):
        considered = (call.coverage is None) or (call.coverage >= min_coverage)
        if call.frac is None:
            methylated = True                      # presence-only file
        else:
            methylated = call.frac >= min_frac
        hit_idxs = call_hits.get(ci, [])
        motifs_hit = sorted({hits[hi].motif for hi in hit_idxs})
        on_motif = len(hit_idxs) > 0
        mb_key = (call.seqid, "*" if ignore_motif_strand else call.strand, call.pos)
        on_mod_base = bool(modbase.get(mb_key))
        annotations.append(
            CallAnnotation(
                feature_idxs=call_features.get(ci, []),
                hit_idxs=hit_idxs,
                motifs_hit=motifs_hit,
                on_motif=on_motif,
                on_mod_base=on_mod_base,
                considered=considered,
                methylated=methylated,
            )
        )
    return annotations


def assign_motifs_to_features(
    hits: List[MotifHit],
    features: List[Feature],
) -> Dict[int, Dict[str, int]]:
    """feature_idx -> {motif -> occurrence_count} (occurrence located by modified
    base if known, else by motif start)."""
    pts_by_seq: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for hi, h in enumerate(hits):
        rep = h.mod_pos if h.mod_pos is not None else h.start
        pts_by_seq[h.seqid].append((rep, hi))
    feat_by_seq: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for fi, f in enumerate(features):
        feat_by_seq[f.seqid].append((f.start, f.end, fi))
    out: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for seqid, pts in pts_by_seq.items():
        hit_feats = sweep_points_in_intervals(pts, feat_by_seq.get(seqid, []))
        for hi, fidxs in hit_feats.items():
            motif = hits[hi].motif
            for fi in fidxs:
                out[fi][motif] += 1
    return out


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
def w(fh, row):
    fh.writerow(row)


def write_calls_tsv(path, calls, anns, features):
    with open(path, "w", newline="") as f:
        out = csv.writer(f, delimiter="\t")
        out.writerow(["replicon", "pos_0based", "pos_1based", "strand", "mod_code",
                      "coverage", "frac_modified", "considered", "methylated",
                      "n_features", "feature_ids", "feature_types",
                      "on_motif", "on_mod_base", "motifs_hit"])
        for call, a in zip(calls, anns):
            fids = [features[i].fid for i in a.feature_idxs]
            ftys = sorted({features[i].ftype for i in a.feature_idxs})
            out.writerow([
                call.seqid, call.pos, call.pos + 1, call.strand, call.mod_code,
                "" if call.coverage is None else f"{call.coverage:g}",
                "" if call.frac is None else f"{call.frac:.4f}",
                int(a.considered), int(a.methylated),
                len(fids), ",".join(fids), ",".join(ftys),
                int(a.on_motif), int(a.on_mod_base), ",".join(a.motifs_hit),
            ])


def write_per_feature(path, features, calls, anns, hits, feat_motifs, motif_list):
    """
    Aggregate per feature. A methylation call contributes to a feature's
    'on motif' counts only when ALL of the following lie inside the feature
    [start, end]:
        * the methylation position (guaranteed: the call is assigned to the
          feature by genomic containment),
        * the modified/related base of the motif occurrence the call sits on
          (or, when no offset is provided, the motif occurrence's start),
    and the call is on that motif occurrence (strand-aware).
    This makes the three-way (methylation / motif base / feature) overlap explicit.
    """
    agg = {fi: dict(n_calls=0, n_considered=0, n_meth=0,
                    n_calls_on_motif=0, n_meth_on_motif=0, n_meth_on_modbase=0)
           for fi in range(len(features))}
    for call, a in zip(calls, anns):
        for fi in a.feature_idxs:
            feat = features[fi]
            d = agg[fi]
            d["n_calls"] += 1
            # is the call on a motif occurrence whose related base sits in THIS feature?
            on_motif_here = False
            on_modbase_here = False
            for hi in a.hit_idxs:
                h = hits[hi]
                rep = h.mod_pos if h.mod_pos is not None else h.start
                if feat.start <= rep <= feat.end:
                    on_motif_here = True
                    if h.mod_pos is not None and h.mod_pos == call.pos:
                        on_modbase_here = True
            if a.considered:
                d["n_considered"] += 1
                if a.methylated:
                    d["n_meth"] += 1
                    if on_motif_here:
                        d["n_meth_on_motif"] += 1
                    if on_modbase_here:
                        d["n_meth_on_modbase"] += 1
                if on_motif_here:
                    d["n_calls_on_motif"] += 1
    with open(path, "w", newline="") as f:
        out = csv.writer(f, delimiter="\t")
        out.writerow(["feature_id", "replicon", "type", "strand",
                      "start_1based", "end_1based", "length",
                      "n_motif_occurrences", "motifs_present", "has_motif",
                      "n_calls", "n_considered", "n_methylated",
                      "n_calls_on_motif", "n_methylated_on_motif",
                      "n_methylated_on_modbase", "n_methylated_off_motif",
                      "frac_methylated_on_motif", "is_methylated_feature",
                      "methylation_on_motif"])
        for fi, feat in enumerate(features):
            d = agg[fi]
            fm = feat_motifs.get(fi, {})
            n_occ = sum(fm.values())
            present = sorted(m for m in motif_list if fm.get(m, 0) > 0)
            n_meth = d["n_meth"]
            n_on = d["n_meth_on_motif"]
            n_off = n_meth - n_on
            frac_on = (n_on / n_meth) if n_meth else 0.0
            out.writerow([
                feat.fid, feat.seqid, feat.ftype, feat.strand,
                feat.start + 1, feat.end + 1, feat.end - feat.start + 1,
                n_occ, ",".join(present), int(n_occ > 0),
                d["n_calls"], d["n_considered"], n_meth,
                d["n_calls_on_motif"], n_on, d["n_meth_on_modbase"], n_off,
                f"{frac_on:.4f}", int(n_meth > 0), int(n_on > 0),
            ])
    return agg


def write_per_feature_motif(path, features, calls, anns, hits, feat_motifs, motif_list):
    # per (feature, motif): occurrences + methylated calls whose motif related
    # base falls inside the feature (same explicit linkage as write_per_feature)
    meth_on_motif = defaultdict(lambda: defaultdict(int))  # fi -> motif -> count
    for call, a in zip(calls, anns):
        if not (a.considered and a.methylated and a.hit_idxs):
            continue
        for fi in a.feature_idxs:
            feat = features[fi]
            for hi in a.hit_idxs:
                h = hits[hi]
                rep = h.mod_pos if h.mod_pos is not None else h.start
                if feat.start <= rep <= feat.end:
                    meth_on_motif[fi][h.motif] += 1
    with open(path, "w", newline="") as f:
        out = csv.writer(f, delimiter="\t")
        out.writerow(["feature_id", "replicon", "type", "motif",
                      "n_motif_occurrences", "n_methylated_calls_on_motif"])
        for fi, feat in enumerate(features):
            fm = feat_motifs.get(fi, {})
            for m in motif_list:
                occ = fm.get(m, 0)
                mm = meth_on_motif.get(fi, {}).get(m, 0)
                if occ == 0 and mm == 0:
                    continue
                out.writerow([feat.fid, feat.seqid, feat.ftype, m, occ, mm])


def write_motif_summary(path, hits, calls, anns, motif_list):
    occ = defaultdict(int)
    for h in hits:
        occ[h.motif] += 1
    meth_on = defaultdict(int)
    calls_on = defaultdict(int)
    modbase_meth = defaultdict(int)
    for call, a in zip(calls, anns):
        if not a.considered:
            continue
        for m in a.motifs_hit:
            calls_on[m] += 1
            if a.methylated:
                meth_on[m] += 1
        if a.on_mod_base and a.methylated:
            for m in a.motifs_hit:
                modbase_meth[m] += 1
    with open(path, "w", newline="") as f:
        out = csv.writer(f, delimiter="\t")
        out.writerow(["motif", "n_occurrences_both_strands",
                      "n_assayed_calls_on_motif", "n_methylated_calls_on_motif",
                      "n_methylated_on_modbase",
                      "frac_assayed_on_motif_methylated"])
        for m in motif_list:
            con = calls_on[m]
            frac = (meth_on[m] / con) if con else 0.0
            out.writerow([m, occ[m], con, meth_on[m], modbase_meth[m], f"{frac:.4f}"])


# --------------------------------------------------------------------------- #
# Enrichment
# --------------------------------------------------------------------------- #
def compute_enrichment(calls, anns, features, feat_agg, feat_motifs,
                       motif_list, offsets_known, feature_meth_min_calls):
    """
    Build the two 2x2 contingency tests per motif and return a list of dicts.
    Odds ratio is computed directly (Haldane-corrected) so it is always
    available; the exact p-value uses scipy when present, else NaN. BH q-values
    are added within each level.
    """
    try:
        from scipy.stats import fisher_exact
        have_scipy = True
    except ImportError:
        have_scipy = False
        sys.stderr.write("[warn] scipy not available; enrichment p-values = NaN.\n")

    def odds_ratio(a, b, c, d):
        # Haldane-Anscombe 0.5 correction to keep it finite/interpretable
        return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

    rows = []
    # ---- Site level: among assayed calls, on-motif-site vs methylated ----
    for m in motif_list:
        a = b = c = d = 0
        for ann in anns:
            if not ann.considered:
                continue
            if offsets_known:
                on_site = (m in ann.motifs_hit) and ann.on_mod_base
            else:
                on_site = (m in ann.motifs_hit)
            meth = ann.methylated
            if on_site and meth:
                a += 1
            elif on_site and not meth:
                b += 1
            elif (not on_site) and meth:
                c += 1
            else:
                d += 1
        p = float("nan")
        if have_scipy and (a + b) > 0:
            try:
                _, p = fisher_exact([[a, b], [c, d]])
            except ValueError:
                p = float("nan")
        rows.append(dict(level="site", motif=m, a=a, b=b, c=c, d=d,
                         odds=odds_ratio(a, b, c, d), p=p))

    # ---- Feature level: among features, has-motif vs feature-methylated ----
    for m in motif_list:
        a = b = c = d = 0
        for fi in range(len(features)):
            has_m = feat_motifs.get(fi, {}).get(m, 0) > 0
            meth_feat = feat_agg[fi]["n_meth"] >= feature_meth_min_calls
            if has_m and meth_feat:
                a += 1
            elif has_m and not meth_feat:
                b += 1
            elif (not has_m) and meth_feat:
                c += 1
            else:
                d += 1
        p = float("nan")
        if have_scipy:
            try:
                _, p = fisher_exact([[a, b], [c, d]])
            except ValueError:
                p = float("nan")
        rows.append(dict(level="feature", motif=m, a=a, b=b, c=c, d=d,
                         odds=odds_ratio(a, b, c, d), p=p))

    # BH within each level
    by_level = defaultdict(list)
    for idx, r in enumerate(rows):
        by_level[r["level"]].append(idx)
    for _level, idxs in by_level.items():
        ps = [rows[i]["p"] for i in idxs]
        qs = benjamini_hochberg([p if p == p else 1.0 for p in ps])
        for i, q in zip(idxs, qs):
            rows[i]["q"] = q
    return rows


def write_enrichment(path, rows):
    with open(path, "w", newline="") as f:
        out = csv.writer(f, delimiter="\t")
        out.writerow(["level", "motif",
                      "a_onmotif_meth", "b_onmotif_unmeth",
                      "c_offmotif_meth", "d_offmotif_unmeth",
                      "odds_ratio", "p_value", "q_value_BH"])
        for r in rows:
            out.writerow([r["level"], r["motif"], r["a"], r["b"], r["c"], r["d"],
                          f"{r['odds']:.4g}", f"{r['p']:.4g}",
                          f"{r.get('q', float('nan')):.4g}"])


# --------------------------------------------------------------------------- #
# Auxiliary data dumps (decoupled re-plotting / external tools)
# --------------------------------------------------------------------------- #
def write_replicons_tsv(path, genome):
    with open(path, "w", newline="") as f:
        out = csv.writer(f, delimiter="\t")
        out.writerow(["replicon", "length"])
        for k in sorted(genome, key=lambda x: -len(genome[x])):
            out.writerow([k, len(genome[k])])


def write_motif_hits_tsv(path, hits):
    with open(path, "w", newline="") as f:
        out = csv.writer(f, delimiter="\t")
        out.writerow(["replicon", "start_0based", "end_0based", "strand",
                      "motif", "mod_pos_0based"])
        for h in hits:
            out.writerow([h.seqid, h.start, h.end, h.strand, h.motif,
                          "" if h.mod_pos is None else h.mod_pos])


def read_fasta_labels(path, strip_version):
    """Build short, human-readable replicon labels from FASTA headers
    (e.g. 'chromosome', 'pSymA') for nicer plot annotation."""
    labels = {}
    if not path:
        return labels
    try:
        with smart_open(path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                sid = header.split()[0]
                if strip_version:
                    sid = _strip_version(sid)
                desc = header.lower()
                lab = sid
                if "chromosome" in desc:
                    lab = "chromosome"
                elif "plasmid" in desc:
                    after = desc.split("plasmid", 1)[1].strip()
                    tok = after.split(",")[0].split()[0] if after else ""
                    lab = tok if tok else "plasmid"
                labels[sid] = lab
    except OSError:
        pass
    return labels


# --------------------------------------------------------------------------- #
# Visualisation
# --------------------------------------------------------------------------- #
def make_plots(outdir, prefix, genome, features, calls, anns, hits, motif_list,
               feat_agg, feat_motifs, enrich_rows, labels, plot_format, min_frac):
    """Generate insightful figures. Returns the list of written file paths.

    Figures:
      * <prefix>.circular.<motif>.<fmt>  -- motif occurrences over the genome
        (outer ring), with methylated-on-motif density (inner ring) and the
        replicon ideogram in between. One per motif.
      * <prefix>.motif_methylation_summary.<fmt>  -- methylation rate per motif,
        site-level enrichment (odds ratio), and on- vs off-motif fraction.
      * <prefix>.feature_motif_association.<fmt>  -- P(feature methylated) given
        motif presence vs absence, per motif (the 'driven by motif' view).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
        import numpy as np
    except ImportError:
        sys.stderr.write("[warn] matplotlib/numpy not available; skipping --plots. "
                         "Install with: pip install matplotlib numpy\n")
        return []

    written = []
    rep_order = sorted(genome, key=lambda x: -len(genome[x]))
    rep_len = {k: len(genome[k]) for k in rep_order}
    total = sum(rep_len.values()) or 1
    rep_color = {k: plt.cm.tab10(i % 10) for i, k in enumerate(rep_order)}
    lab = lambda k: labels.get(k, k)  # noqa: E731

    # angular layout (start at top, clockwise), proportional + gaps
    gap = 2 * math.pi * 0.015
    usable = 2 * math.pi - len(rep_order) * gap
    layout = {}
    cur = math.pi / 2  # 12 o'clock
    for k in rep_order:
        span = usable * rep_len[k] / total
        layout[k] = (cur, cur - span)  # clockwise -> decreasing angle
        cur -= span + gap
    binbp = max(3000, total // 800)

    def hist_theta(positions_by_rep):
        """Return per-replicon (theta_centers, widths, counts) and global max."""
        out = {}
        gmax = 0
        for k in rep_order:
            L = rep_len[k]
            nb = max(2, int(math.ceil(L / binbp)))
            counts, edges = np.histogram(positions_by_rep.get(k, []),
                                         bins=nb, range=(0, L))
            t0, t1 = layout[k]
            theta_edges = t0 + (edges / L) * (t1 - t0)
            centers = (theta_edges[:-1] + theta_edges[1:]) / 2
            widths = np.abs(np.diff(theta_edges))
            out[k] = (centers, widths, counts)
            if counts.size:
                gmax = max(gmax, counts.max())
        return out, max(gmax, 1)

    # ---------- circular plot per motif ----------
    for m in motif_list:
        occ = defaultdict(list)
        for h in hits:
            if h.motif == m:
                occ[h.seqid].append(h.mod_pos if h.mod_pos is not None else h.start)
        meth = defaultdict(list)
        for call, a in zip(calls, anns):
            if a.considered and a.methylated and (m in a.motifs_hit):
                meth[call.seqid].append(call.pos)

        occ_h, occ_max = hist_theta(occ)
        meth_h, meth_max = hist_theta(meth)

        fig = plt.figure(figsize=(8.2, 8.2))
        ax = fig.add_subplot(111, projection="polar")
        ax.set_theta_zero_location("E")
        ax.set_theta_direction(1)
        ax.axis("off")

        R_IDEO, W_IDEO = 1.00, 0.10
        R_OCC, H_OCC = R_IDEO + W_IDEO + 0.03, 0.62      # outer: occurrences
        R_METH, H_METH = 0.42, 0.50                       # inner: methylated
        for k in rep_order:
            t0, t1 = layout[k]
            tc = (t0 + t1) / 2
            ax.bar(tc, W_IDEO, width=abs(t1 - t0), bottom=R_IDEO,
                   color=rep_color[k], edgecolor="white", linewidth=0.6,
                   align="center", zorder=3)
            # occurrences (outward, steelblue)
            c, w, cnt = occ_h[k]
            ax.bar(c, (cnt / occ_max) * H_OCC, width=w, bottom=R_OCC,
                   color="#3a6ea5", linewidth=0, align="center", zorder=2)
            # methylated-on-motif (inward ring, crimson)
            c, w, cnt = meth_h[k]
            ax.bar(c, (cnt / meth_max) * H_METH, width=w, bottom=R_METH,
                   color="#c0392b", linewidth=0, align="center", zorder=2)
            # replicon label
            deg = math.degrees(tc) % 360
            ha = "left" if (deg < 90 or deg > 270) else "right"
            ax.text(tc, R_OCC + H_OCC + 0.12, f"{lab(k)}\n{rep_len[k]/1e6:.2f} Mb",
                    rotation=0, ha=ha, va="center", fontsize=9,
                    color=rep_color[k], fontweight="bold")

        ax.set_ylim(0, R_OCC + H_OCC + 0.30)
        n_occ_tot = sum(len(v) for v in occ.values())
        n_meth_tot = sum(len(v) for v in meth.values())
        ax.set_title(f"{m}   |   {n_occ_tot} occurrences   |   "
                     f"{n_meth_tot} methylated on motif",
                     fontsize=13, fontweight="bold", pad=18)
        ax.text(0, 0, m, ha="center", va="center", fontsize=11,
                fontweight="bold", color="#333333")
        legend = [Patch(facecolor="#3a6ea5", label="motif occurrences (outer)"),
                  Patch(facecolor="#c0392b", label="methylated on motif (inner)")]
        ax.legend(handles=legend, loc="lower center",
                  bbox_to_anchor=(0.5, -0.06), ncol=2, frameon=False, fontsize=9)
        path = os.path.join(outdir, f"{prefix}.circular.{m}.{plot_format}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    # ---------- summary figure ----------
    # methylation rate per motif (assayed on-motif sites that are methylated)
    on_assayed = {m: 0 for m in motif_list}
    on_meth = {m: 0 for m in motif_list}
    off_frac, on_frac = [], {m: [] for m in motif_list}
    for call, a in zip(calls, anns):
        if not a.considered or call.frac is None:
            continue
        if a.motifs_hit:
            for m in a.motifs_hit:
                on_assayed[m] += 1
                if a.methylated:
                    on_meth[m] += 1
                on_frac[m].append(call.frac)
        else:
            off_frac.append(call.frac)
    rate = [(on_meth[m] / on_assayed[m]) if on_assayed[m] else 0.0
            for m in motif_list]

    site = {r["motif"]: r for r in enrich_rows if r["level"] == "site"}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    x = np.arange(len(motif_list))
    axes[0].bar(x, rate, color="#c0392b")
    axes[0].set_xticks(x); axes[0].set_xticklabels(motif_list, rotation=30, ha="right")
    axes[0].set_ylabel("fraction methylated"); axes[0].set_ylim(0, 1)
    axes[0].set_title("Methylation rate of assayed on-motif sites")
    for xi, ra, m in zip(x, rate, motif_list):
        axes[0].text(xi, ra + 0.02, f"{ra:.0%}\n(n={on_assayed[m]})",
                     ha="center", va="bottom", fontsize=8)

    ors = [site[m]["odds"] if m in site else float("nan") for m in motif_list]
    qs = [site[m].get("q", float("nan")) if m in site else float("nan")
          for m in motif_list]
    colors = ["#3a6ea5" if o >= 1 else "#7f8c8d" for o in ors]
    axes[1].bar(x, ors, color=colors)
    axes[1].axhline(1.0, color="black", lw=0.8, ls="--")
    axes[1].set_yscale("log")
    axes[1].set_xticks(x); axes[1].set_xticklabels(motif_list, rotation=30, ha="right")
    axes[1].set_ylabel("odds ratio (log)")
    axes[1].set_title("Site-level enrichment\n(methylated ~ on motif)")
    for xi, o, q in zip(x, ors, qs):
        star = "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "ns"
        axes[1].text(xi, o, f" {star}", ha="center",
                     va="bottom" if o >= 1 else "top", fontsize=9)

    data = [on_frac[m] for m in motif_list] + [off_frac]
    box = axes[2].boxplot(data, showfliers=False, patch_artist=True,
                          medianprops=dict(color="black"))
    palette = ["#c0392b"] * len(motif_list) + ["#bdc3c7"]
    for patch, col in zip(box["boxes"], palette):
        patch.set_facecolor(col)
    axes[2].axhline(min_frac, color="#888", lw=0.8, ls=":")
    axes[2].set_xticklabels(motif_list + ["off-motif"], rotation=30, ha="right")
    axes[2].set_ylabel("fraction modified per call")
    axes[2].set_title("On-motif vs background methylation")
    fig.tight_layout()
    path = os.path.join(outdir, f"{prefix}.motif_methylation_summary.{plot_format}")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    written.append(path)

    # ---------- feature-level association figure ----------
    feat = {r["motif"]: r for r in enrich_rows if r["level"] == "feature"}
    p_with, p_without = [], []
    for m in motif_list:
        r = feat.get(m)
        if r:
            a, b, c, d = r["a"], r["b"], r["c"], r["d"]
            p_with.append(a / (a + b) if (a + b) else 0.0)
            p_without.append(c / (c + d) if (c + d) else 0.0)
        else:
            p_with.append(0.0); p_without.append(0.0)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    width = 0.38
    ax.bar(x - width / 2, p_with, width, label="feature contains motif",
           color="#27ae60")
    ax.bar(x + width / 2, p_without, width, label="feature lacks motif",
           color="#95a5a6")
    ax.set_xticks(x); ax.set_xticklabels(motif_list, rotation=20, ha="right")
    ax.set_ylabel("fraction of features methylated"); ax.set_ylim(0, 1)
    ax.set_title("Is feature methylation associated with motif presence?")
    for xi, m in zip(x, motif_list):
        r = feat.get(m)
        if r:
            q = r.get("q", float("nan"))
            star = "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "ns"
            ax.text(xi, max(p_with[list(motif_list).index(m)],
                            p_without[list(motif_list).index(m)]) + 0.03,
                    f"OR={r['odds']:.1f} {star}", ha="center", fontsize=8)
    ax.legend(frameon=False)
    fig.tight_layout()
    path = os.path.join(outdir, f"{prefix}.feature_motif_association.{plot_format}")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    written.append(path)
    return written


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Cross methylation calls, DNA motifs and GFF3 annotation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--gff", required=True, help="GFF3 annotation (may embed FASTA).")
    p.add_argument("--methylation", required=True,
                   help="Methylation calls (bedMethyl or GFF3).")
    p.add_argument("--motifs", required=True,
                   help="Newline-delimited motifs; optional 2nd column = 0-based "
                        "modified-base offset.")
    p.add_argument("--fasta", default=None,
                   help="Reference FASTA. If omitted, an embedded ##FASTA block in "
                        "the GFF3 is used.")
    p.add_argument("--meth-format", default="auto",
                   choices=["auto", "bedmethyl", "gff"],
                   help="Methylation file format.")
    p.add_argument("--mod-codes", default=None,
                   help="Comma-separated modification codes/types to keep "
                        "(e.g. 'a,m'). Default: keep all.")
    p.add_argument("--feature-types", default=None,
                   help="Comma-separated GFF feature types to keep "
                        "(e.g. 'CDS,gene'). Default: keep all.")
    p.add_argument("--id-attr", default="ID",
                   help="GFF attribute used as feature id (falls back to "
                        "locus_tag/Name, then a synthesized id).")
    p.add_argument("--min-coverage", type=float, default=0.0,
                   help="Minimum valid coverage for a call to be considered.")
    p.add_argument("--min-frac", type=float, default=0.5,
                   help="Fraction-modified threshold to call a position methylated "
                        "(used when the file carries a fraction).")
    p.add_argument("--ignore-motif-strand", action="store_true",
                   help="Match motifs irrespective of strand.")
    p.add_argument("--feature-meth-min-calls", type=int, default=1,
                   help="Min methylated calls for a feature to be 'methylated' "
                        "(feature-level enrichment).")
    p.add_argument("--fisher", action="store_true",
                   help="Run Fisher exact enrichment tests (site & feature level) "
                        "with BH/FDR correction.")
    p.add_argument("--strip-version", action="store_true",
                   help="Strip a trailing .N version suffix from ALL sequence ids "
                        "(GFF, FASTA, methylation) before crossing, to reconcile "
                        "e.g. 'NZ_CP012345.1' vs 'NZ_CP012345'.")
    p.add_argument("--plots", action="store_true",
                   help="Generate figures (circular per-motif genome plots plus "
                        "summary/enrichment views). Requires matplotlib.")
    p.add_argument("--plot-format", default="png", choices=["png", "pdf", "svg"],
                   help="Image format for --plots.")
    p.add_argument("--outdir", default="motif_meth_out", help="Output directory.")
    p.add_argument("--prefix", default="sample", help="Output filename prefix.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


_VERSION_RE = re.compile(r"\.\d+$")


def _strip_version(seqid: str) -> str:
    return _VERSION_RE.sub("", seqid)


def reconcile_seqids(features, genome, calls, strip_version):
    """
    Ensure the three layers share a coordinate namespace. Motif search runs on
    the FASTA seqids, so features and calls must use the SAME seqids or they
    silently fail to cross. If the raw ids are disjoint but stripping a trailing
    .N version suffix makes them agree, this is applied automatically (with a
    loud notice). Returns (genome, strip_applied). Raises SystemExit on a fatal,
    unrecoverable mismatch.
    """
    gff_raw = {f.seqid for f in features}
    fa_raw = set(genome)

    apply_strip = strip_version
    if not apply_strip and gff_raw and not (gff_raw & fa_raw):
        gff_s = {_strip_version(x) for x in gff_raw}
        fa_s = {_strip_version(x) for x in fa_raw}
        if gff_s & fa_s:
            sys.stderr.write(
                "[auto] GFF and FASTA seqids differ only by a version suffix "
                "(e.g. '.1'); stripping it from all layers to reconcile.\n")
            apply_strip = True

    if apply_strip:
        for f in features:
            f.seqid = _strip_version(f.seqid)
        for c in calls:
            c.seqid = _strip_version(c.seqid)
        genome = {_strip_version(k): v for k, v in genome.items()}

    gff_ids = {f.seqid for f in features}
    fa_ids = set(genome)
    me_ids = {c.seqid for c in calls}

    def fmt(s):
        s = sorted(s)
        return ", ".join(s[:4]) + (" ..." if len(s) > 4 else "")

    sys.stderr.write(
        "      seqids -> "
        f"GFF: [{fmt(gff_ids)}] | FASTA: [{fmt(fa_ids)}] | meth: [{fmt(me_ids)}]\n"
    )

    fatal = False
    if not (gff_ids & fa_ids):
        sys.stderr.write(
            "[ERROR] GFF feature seqids do not match any FASTA seqid. Motifs are "
            "located on FASTA coordinates, so NO methylation/motif can be placed "
            "in a feature -> all per-feature counts would be 0. "
            "Use --strip-version, or provide a GFF/FASTA with matching ids.\n")
        fatal = True
    if calls and not (me_ids & fa_ids):
        sys.stderr.write(
            "[ERROR] methylation seqids do not match any FASTA seqid -> no call "
            "can fall on a motif. Use --strip-version or fix the ids.\n")
        fatal = True
    if fatal:
        raise SystemExit(3)

    only_meth = sorted(me_ids - fa_ids)
    if only_meth:
        sys.stderr.write(f"[warn] {len(only_meth)} methylation seqid(s) absent from "
                         f"FASTA, e.g. {only_meth[:3]} (those calls are ignored).\n")
    return genome, apply_strip


def main(argv=None) -> int:
    args = build_argparser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    feature_types = set(args.feature_types.split(",")) if args.feature_types else None
    keep_codes = set(args.mod_codes.split(",")) if args.mod_codes else None

    sys.stderr.write("[1/6] Parsing GFF3 ...\n")
    features, embedded = parse_gff3(args.gff, feature_types, args.id_attr)
    sys.stderr.write(f"      {len(features)} features"
                     f"{' (filtered)' if feature_types else ''}.\n")

    sys.stderr.write("[2/6] Loading reference sequence ...\n")
    if args.fasta:
        genome = load_fasta(args.fasta)
    elif embedded:
        genome = embedded
        sys.stderr.write("      using FASTA embedded in GFF3.\n")
    else:
        sys.stderr.write("ERROR: no sequence available. Provide --fasta or embed "
                         "a ##FASTA block in the GFF3 (motif search needs the "
                         "genome sequence).\n")
        return 2
    sys.stderr.write(f"      {len(genome)} replicon(s), "
                     f"{sum(len(s) for s in genome.values())} bp.\n")

    sys.stderr.write("[3/6] Parsing motifs and methylation, reconciling seqids ...\n")
    motifs = parse_motifs(args.motifs)
    offsets_known = all(off is not None for _m, off in motifs)
    motif_list = [m for m, _ in motifs]
    calls = parse_methylation(args.methylation, args.meth_format, keep_codes)
    sys.stderr.write(f"      {len(motifs)} motif(s), {len(calls)} methylation "
                     f"call(s). modified-base offsets "
                     f"{'provided' if offsets_known else 'NOT fully provided'}.\n")

    # CRITICAL: features, FASTA and methylation must share a seqid namespace,
    # otherwise the cross silently yields all-zero per-feature counts.
    genome, _strip_applied = reconcile_seqids(features, genome, calls,
                                              args.strip_version)

    sys.stderr.write("[4/6] Searching motifs on the genome ...\n")
    hits = find_motif_hits(genome, motifs)
    sys.stderr.write(f"      {len(hits)} motif occurrence(s) (both strands).\n")

    sys.stderr.write("[5/6] Crossing layers ...\n")
    anns = annotate_calls(calls, features, hits,
                          args.min_coverage, args.min_frac,
                          args.ignore_motif_strand, offsets_known)
    feat_motifs = assign_motifs_to_features(hits, features)

    sys.stderr.write("[6/6] Writing outputs ...\n")
    pfx = os.path.join(args.outdir, args.prefix)
    write_calls_tsv(f"{pfx}.methylation_calls.tsv", calls, anns, features)
    feat_agg = write_per_feature(f"{pfx}.per_feature.tsv", features, calls, anns,
                                 hits, feat_motifs, motif_list)
    write_per_feature_motif(f"{pfx}.per_feature_motif.tsv", features, calls, anns,
                            hits, feat_motifs, motif_list)
    write_motif_summary(f"{pfx}.motif_summary.tsv", hits, calls, anns, motif_list)
    write_replicons_tsv(f"{pfx}.replicons.tsv", genome)
    write_motif_hits_tsv(f"{pfx}.motif_hits.tsv", hits)

    # enrichment is needed for the figures too, so compute when either is requested
    enrich_rows = []
    if args.fisher or args.plots:
        enrich_rows = compute_enrichment(calls, anns, features, feat_agg,
                                         feat_motifs, motif_list, offsets_known,
                                         args.feature_meth_min_calls)
    if args.fisher:
        write_enrichment(f"{pfx}.enrichment.tsv", enrich_rows)

    if args.plots:
        sys.stderr.write("      rendering figures ...\n")
        labels = read_fasta_labels(args.fasta, True)
        figs = make_plots(args.outdir, args.prefix, genome, features, calls, anns,
                          hits, motif_list, feat_agg, feat_motifs, enrich_rows,
                          labels, args.plot_format, args.min_frac)
        for fp in figs:
            sys.stderr.write(f"        {fp}\n")

    # ---- console summary ----
    n_considered = sum(1 for a in anns if a.considered)
    n_meth = sum(1 for a in anns if a.considered and a.methylated)
    n_meth_on = sum(1 for a in anns if a.considered and a.methylated and a.on_motif)
    n_feat_meth = sum(1 for fi in range(len(features))
                      if feat_agg[fi]["n_meth"] >= args.feature_meth_min_calls)
    n_feat_motif = sum(1 for fi in range(len(features))
                       if sum(feat_motifs.get(fi, {}).values()) > 0)
    pct = (100 * n_meth_on / n_meth) if n_meth else 0.0
    sys.stderr.write(
        "\n=== Summary ===\n"
        f"  calls considered (cov>= {args.min_coverage:g}): {n_considered}\n"
        f"  methylated (frac>= {args.min_frac:g} or presence): {n_meth}\n"
        f"  methylated AND on a motif: {n_meth_on} ({pct:.1f}% of methylated)\n"
        f"  features with >=1 motif occurrence: {n_feat_motif}/{len(features)}\n"
        f"  features methylated (>= {args.feature_meth_min_calls} calls): "
        f"{n_feat_meth}/{len(features)}\n"
        f"  outputs written to: {args.outdir}/\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
