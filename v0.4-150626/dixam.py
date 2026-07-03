#!/usr/bin/env python3
"""
dixam.py
---------------------------
Differential version of xam_2.py.

Cross-references a Prokka GFF3 annotation with methylation basecalling from
TWO conditions of the SAME sample — a control and a treatment — and reports,
per annotated feature and per replicon, how methylation is distributed in
each condition and how it differs between them.

Inputs
------
    -g / --gff                 Prokka GFF3 annotation (shared by both conditions).
    Control methylation (exactly one):
        -b / --bed             Filtered 6-column BED (control).
        -m / --methyl-gff      Basecalling GFF3 (control).
    Treatment methylation (exactly one):
        -B / --bed2            Filtered 6-column BED (treatment).
        -M / --methyl-gff2     Basecalling GFF3 (treatment).

For every annotated feature × methylation type the script reports, for each
condition, the number of methylated positions inside the feature (raw_count),
and derives the differential metrics described below.

Per-condition output ("listed singularly per BED file")
-------------------------------------------------------
    <stem>_<control_label>_<MOD>.tsv
    <stem>_<treatment_label>_<MOD>.tsv
        Same per-feature listing as xam_2, one set per condition.

Crossed differential output
---------------------------
    <stem>_differential_<MOD>.tsv   (or a single file with --no-split)
        One row per feature × mod where either condition has >0 calls, with:
            <ctrl>_raw / <trt>_raw      raw methylated-position counts
            <ctrl>_cpm / <trt>_cpm      depth-normalised counts (calls per
                                        million total calls in that condition)
            <ctrl>_norm / <trt>_norm    raw_count / feature_length_bp
            delta_raw                   trt_raw - ctrl_raw
            delta_cpm                   trt_cpm - ctrl_cpm
            log2FC_cpm                  log2((trt_cpm+pc)/(ctrl_cpm+pc))
            pvalue / qvalue             per-feature Fisher exact test on the
                                        2x2 table of in-feature vs remaining
                                        calls of that mod (control vs
                                        treatment), BH-adjusted across features
                                        (only if scipy is available; disable
                                        with --no-test)

Normalisation note
------------------
The 6-column filtered BED only carries methylated-base POSITIONS, not
coverage, so this tool compares the *distribution* of methylation calls, not
methylation fraction. Conditions are made comparable by library-size scaling
(CPM = count / total_calls_in_condition x 1e6). For coverage-aware, fraction-
based differential methylation (DMR calling) a dedicated tool (methylKit, DSS,
modkit dmr) operating on full bedMethyl is required; treat the qvalues here as
a screen for differentially-distributed features, not a calibrated DMR result.

Figures
-------
    <stem>_along_replicons_overlay.<fmt>
        Methylation distribution along each replicon (CPM per positional bin),
        control vs treatment overlaid, one row per replicon and one column per
        methylation type.
    <stem>_along_replicons_diff.<fmt>
        The same axis showing (treatment - control) CPM per bin as a diverging
        track: positive (gain in treatment) vs negative (loss).
    <stem>_per_replicon_feature.<fmt>
        Only when more than one feature type carries methylation: control vs
        treatment CPM per feature type, grouped bars, one subplot per
        replicon x methylation type.

Usage
-----
    python3 dixam.py -g annotation.gff \\
        -b control.bed -B treatment.bed -o results/

    python3 dixam.py -g annotation.gff \\
        -m control.basecall.gff -B treatment.bed -o results/ \\
        --labels WT secDF_kd --features CDS tRNA rRNA
"""

import argparse
import logging
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MOD_CODE_TO_LABEL: dict[str, str] = {
    "m": "5mC",
    "h": "5hmC",
    "a": "6mA",
    "f": "5fC",
    "c": "5caC",
    "g": "7mG",
    "e": "3mC",
    "b": "5hmU",
}

MOD_LABEL_TO_DESC: dict[str, str] = {
    "5mC":  "5-methylcytosine",
    "m5C":  "5-methylcytosine",
    "m4C":  "4-methylcytosine",
    "5hmC": "5-hydroxymethylcytosine",
    "6mA":  "6-methyladenine",
    "m6A":  "6-methyladenine",
    "5fC":  "5-formylcytosine",
    "5caC": "5-carboxylcytosine",
    "7mG":  "N7-methylguanine",
    "3mC":  "N3-methylcytosine",
    "5hmU": "5-hydroxymethyluracil",
}

MOD_LABEL_TO_CODE: dict[str, str] = {v: k for k, v in MOD_CODE_TO_LABEL.items()}

METHYL_GFF_TYPE_SYNONYMS: dict[str, str] = {
    "cpg":           "m",
    "5mc":           "m",
    "m5C":           "m",
    "5-mc":          "m",
    "methylation":   "m",
    "modified_base": "m",
    "6ma":           "a",
    "m6A":           "a",
    "6-ma":          "a",
    "dam":           "a",
    "gatc":          "a",
    "dcm":           "m",
    "ccwgg":         "m",
    "5hmc":          "h",
    "5-hmc":         "h",
}

# GFF3 column indices (0-based)
GFF_SEQID  = 0
GFF_TYPE   = 2
GFF_START  = 3
GFF_END    = 4
GFF_STRAND = 6
GFF_ATTRS  = 8

# BED column indices (0-based)
BED_CHROM  = 0
BED_START  = 1
BED_MOD    = 3
BED_STRAND = 5

GFF_META_TYPES = {"region", "sequence_feature", "sequence_alteration"}

# Fixed condition colours used across all figures
CONDITION_PALETTE = ["#2c7fb8", "#de2d26"]   # control = blue, treatment = red

# Per-condition single-listing header (identical to xam/xam_2)
OUTPUT_HEADER = "\t".join([
    "replicon",
    "feature_type",
    "feat_start",
    "feat_end",
    "strand",
    "feat_length_bp",
    "locus_tag",
    "gene",
    "product",
    "mod_code",
    "mod_label",
    "mod_description",
    "raw_count",
    "norm_count",
])


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GffFeature:
    """One annotated feature from a Prokka GFF3 file."""
    seqid:      str
    feat_type:  str
    start:      int    # 0-based
    end:        int    # 0-based exclusive
    strand:     str
    locus_tag:  str = ""
    gene:       str = ""
    product:    str = ""

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class MethylSite:
    """One methylation site (from BED or basecalling GFF)."""
    chrom:  str
    start:  int
    mod:    str
    strand: str


@dataclass
class OutputRow:
    """One per-condition result row: one feature × one mod type."""
    feature:   GffFeature
    mod_code:  str
    raw_count: int

    @property
    def norm_count(self) -> float:
        if self.feature.length == 0:
            return 0.0
        return self.raw_count / self.feature.length

    def to_tsv(self) -> str:
        f     = self.feature
        label = MOD_CODE_TO_LABEL.get(self.mod_code, self.mod_code)
        desc  = MOD_LABEL_TO_DESC.get(label, "unknown modification")
        return "\t".join([
            f.seqid, f.feat_type, str(f.start), str(f.end), f.strand,
            str(f.length), f.locus_tag, f.gene, f.product,
            self.mod_code, label, desc,
            str(self.raw_count), f"{self.norm_count:.8f}",
        ])


@dataclass
class DiffRow:
    """One crossed differential result row: feature × mod, both conditions."""
    feature:   GffFeature
    mod_code:  str
    ctrl_raw:  int = 0
    trt_raw:   int = 0
    ctrl_cpm:  float = 0.0
    trt_cpm:   float = 0.0
    log2fc:    float = 0.0
    pvalue:    float = float("nan")
    qvalue:    float = float("nan")

    @property
    def ctrl_norm(self) -> float:
        return self.ctrl_raw / self.feature.length if self.feature.length else 0.0

    @property
    def trt_norm(self) -> float:
        return self.trt_raw / self.feature.length if self.feature.length else 0.0

    def to_tsv(self, do_test: bool) -> str:
        f     = self.feature
        label = MOD_CODE_TO_LABEL.get(self.mod_code, self.mod_code)
        desc  = MOD_LABEL_TO_DESC.get(label, "unknown modification")
        vals  = [
            f.seqid, f.feat_type, str(f.start), str(f.end), f.strand,
            str(f.length), f.locus_tag, f.gene, f.product,
            self.mod_code, label, desc,
            str(self.ctrl_raw), str(self.trt_raw),
            f"{self.ctrl_cpm:.4f}", f"{self.trt_cpm:.4f}",
            f"{self.ctrl_norm:.8f}", f"{self.trt_norm:.8f}",
            str(self.trt_raw - self.ctrl_raw),
            f"{self.trt_cpm - self.ctrl_cpm:.4f}",
            f"{self.log2fc:.4f}",
        ]
        if do_test:
            vals.append(f"{self.pvalue:.3e}" if self.pvalue == self.pvalue else "NA")
            vals.append(f"{self.qvalue:.3e}" if self.qvalue == self.qvalue else "NA")
        return "\t".join(vals)


def diff_header(labels: tuple[str, str], do_test: bool) -> str:
    c, t = labels
    cols = [
        "replicon", "feature_type", "feat_start", "feat_end", "strand",
        "feat_length_bp", "locus_tag", "gene", "product",
        "mod_code", "mod_label", "mod_description",
        f"{c}_raw", f"{t}_raw",
        f"{c}_cpm", f"{t}_cpm",
        f"{c}_norm", f"{t}_norm",
        "delta_raw", "delta_cpm", "log2FC_cpm",
    ]
    if do_test:
        cols += ["pvalue", "qvalue"]
    return "\t".join(cols)


# ---------------------------------------------------------------------------
# Shared GFF3 parsing
# ---------------------------------------------------------------------------

def _parse_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            key, _, val = part.partition("=")
            attrs[key.strip()] = val.strip()
    return attrs


def _split_gff_line(raw: str) -> list[str] | None:
    cols = raw.split("\t")
    if len(cols) < 9:
        cols = raw.split()
    if len(cols) < 9:
        return None
    return cols


# ---------------------------------------------------------------------------
# Annotation GFF loader  (Prokka)
# ---------------------------------------------------------------------------

def load_annotation_gff(
    gff_path: Path,
    feature_types: set[str] | None = None,
) -> list[GffFeature]:
    features:     list[GffFeature] = []
    skipped_type: int = 0

    with open(gff_path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if raw.startswith("##FASTA"):
                break
            if raw.startswith("#") or not raw:
                continue

            cols = _split_gff_line(raw)
            if cols is None:
                logger.warning(
                    "Annotation GFF line %d: expected 9 columns — skipped.", line_no
                )
                continue

            feat_type = cols[GFF_TYPE]
            if feature_types is not None:
                if feat_type not in feature_types:
                    skipped_type += 1
                    continue
            else:
                if feat_type in GFF_META_TYPES:
                    skipped_type += 1
                    continue

            try:
                start_0 = int(cols[GFF_START]) - 1
                end_0   = int(cols[GFF_END])
            except ValueError:
                logger.warning(
                    "Annotation GFF line %d: cannot parse coordinates — skipped.",
                    line_no,
                )
                continue

            attrs     = _parse_attributes(cols[GFF_ATTRS])
            locus_tag = attrs.get("locus_tag", attrs.get("ID", ""))
            gene      = attrs.get("gene", "")
            product   = attrs.get("product", "")

            features.append(GffFeature(
                seqid     = cols[GFF_SEQID],
                feat_type = feat_type,
                start     = start_0,
                end       = end_0,
                strand    = cols[GFF_STRAND],
                locus_tag = locus_tag,
                gene      = gene,
                product   = product,
            ))

    logger.info(
        "Annotation GFF: %d feature(s) loaded  |  %d skipped by type filter.",
        len(features), skipped_type,
    )
    return features


# ---------------------------------------------------------------------------
# Methylation loaders
# ---------------------------------------------------------------------------

def load_bed(bed_path: Path) -> list[MethylSite]:
    sites:   list[MethylSite] = []
    skipped: int = 0

    with open(bed_path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if not raw or raw.startswith("#"):
                continue

            cols = raw.split("\t")
            if len(cols) < 6:
                logger.warning(
                    "BED line %d: expected >= 6 columns, got %d — skipped.",
                    line_no, len(cols),
                )
                skipped += 1
                continue

            try:
                start = int(cols[BED_START])
            except ValueError:
                logger.warning(
                    "BED line %d: cannot parse start coordinate — skipped.", line_no
                )
                skipped += 1
                continue

            sites.append(MethylSite(
                chrom  = cols[BED_CHROM],
                start  = start,
                mod    = cols[BED_MOD],
                strand = cols[BED_STRAND],
            ))

    logger.info(
        "BED: %d methylation site(s) loaded  |  %d skipped (malformed).",
        len(sites), skipped,
    )
    return sites


def _resolve_mod_code(feat_type: str, attrs: dict[str, str]) -> str:
    if feat_type in MOD_CODE_TO_LABEL:
        return feat_type
    if feat_type in MOD_LABEL_TO_CODE:
        return MOD_LABEL_TO_CODE[feat_type]

    lower = feat_type.lower()
    if lower in METHYL_GFF_TYPE_SYNONYMS:
        return METHYL_GFF_TYPE_SYNONYMS[lower]

    for attr_key in ("mod", "modification", "Name", "type"):
        val = attrs.get(attr_key, "")
        if val in MOD_CODE_TO_LABEL:
            return val
        if val in MOD_LABEL_TO_CODE:
            return MOD_LABEL_TO_CODE[val]
        val_lower = val.lower()
        if val_lower in METHYL_GFF_TYPE_SYNONYMS:
            return METHYL_GFF_TYPE_SYNONYMS[val_lower]

    return feat_type


def load_methyl_gff(gff_path: Path) -> list[MethylSite]:
    sites:         list[MethylSite] = []
    skipped_bad:   int = 0
    skipped_meta:  int = 0
    unknown_types: set[str] = set()

    with open(gff_path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if raw.startswith("##FASTA"):
                break
            if raw.startswith("#") or not raw:
                continue

            cols = _split_gff_line(raw)
            if cols is None:
                logger.warning(
                    "Methylation GFF line %d: expected 9 columns — skipped.", line_no
                )
                skipped_bad += 1
                continue

            feat_type = cols[GFF_TYPE]
            if feat_type in GFF_META_TYPES:
                skipped_meta += 1
                continue

            attrs    = _parse_attributes(cols[GFF_ATTRS])
            mod_code = _resolve_mod_code(feat_type, attrs)

            if mod_code not in MOD_CODE_TO_LABEL and mod_code not in unknown_types:
                unknown_types.add(mod_code)
                logger.warning(
                    "Methylation GFF: unrecognised type '%s' — kept as-is.", mod_code
                )

            try:
                start_0 = int(cols[GFF_START]) - 1
            except ValueError:
                logger.warning(
                    "Methylation GFF line %d: cannot parse start — skipped.", line_no
                )
                skipped_bad += 1
                continue

            sites.append(MethylSite(
                chrom  = cols[GFF_SEQID],
                start  = start_0,
                mod    = mod_code,
                strand = cols[GFF_STRAND],
            ))

    logger.info(
        "Methylation GFF: %d site(s) loaded  |  "
        "%d skipped (malformed)  |  %d skipped (meta-features).",
        len(sites), skipped_bad, skipped_meta,
    )
    if unknown_types:
        logger.warning(
            "Unrecognised modification type(s) kept as raw strings: %s",
            ", ".join(sorted(unknown_types)),
        )
    return sites


def loader_for(bed_path: str | None, gff_path: str | None):
    """Return (resolved_path, source_label, loader_fn) for one condition."""
    if bed_path:
        p = Path(bed_path).resolve()
        return p, f"BED  ({p})", load_bed
    p = Path(gff_path).resolve()
    return p, f"GFF  ({p})", load_methyl_gff


# ---------------------------------------------------------------------------
# Methylation index + overlap counting
# ---------------------------------------------------------------------------

def build_methyl_index(
    sites: list[MethylSite],
) -> dict[tuple[str, str, str], list[int]]:
    index: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    for site in sites:
        index[(site.chrom, site.strand, site.mod)].append(site.start)
    for key in index:
        index[key].sort()
    return index


def count_overlaps(
    feature:       GffFeature,
    index:         dict[tuple[str, str, str], list[int]],
    all_mod_codes: set[str],
    ignore_strand: bool,
) -> dict[str, int]:
    counts: dict[str, int] = {}

    if ignore_strand or feature.strand == ".":
        strands = ["+", "-", "."]
    else:
        strands = [feature.strand]

    for mod in all_mod_codes:
        total = 0
        for strand in strands:
            positions = index.get((feature.seqid, strand, mod))
            if not positions:
                continue

            lo, hi = 0, len(positions)
            while lo < hi:
                mid = (lo + hi) // 2
                if positions[mid] < feature.start:
                    lo = mid + 1
                else:
                    hi = mid

            for pos in positions[lo:]:
                if pos >= feature.end:
                    break
                total += 1

        if total > 0:
            counts[mod] = total

    return counts


def cross_reference(
    features:      list[GffFeature],
    index:         dict[tuple[str, str, str], list[int]],
    all_mod_codes: set[str],
    ignore_strand: bool,
) -> list[OutputRow]:
    """Per-condition listing: one OutputRow per feature × mod (raw_count > 0)."""
    rows: list[OutputRow] = []
    for feat in features:
        counts = count_overlaps(feat, index, all_mod_codes, ignore_strand)
        for mod in sorted(all_mod_codes):
            raw = counts.get(mod, 0)
            if raw > 0:
                rows.append(OutputRow(feature=feat, mod_code=mod, raw_count=raw))
    return rows


def cross_reference_diff(
    features:      list[GffFeature],
    index_ctrl:    dict[tuple[str, str, str], list[int]],
    index_trt:     dict[tuple[str, str, str], list[int]],
    all_mod_codes: set[str],
    ignore_strand: bool,
) -> list[DiffRow]:
    """Crossed listing: one DiffRow per feature × mod where either cond > 0."""
    rows: list[DiffRow] = []
    for feat in features:
        c_counts = count_overlaps(feat, index_ctrl, all_mod_codes, ignore_strand)
        t_counts = count_overlaps(feat, index_trt, all_mod_codes, ignore_strand)
        for mod in sorted(all_mod_codes):
            cr = c_counts.get(mod, 0)
            tr = t_counts.get(mod, 0)
            if cr > 0 or tr > 0:
                rows.append(DiffRow(
                    feature=feat, mod_code=mod, ctrl_raw=cr, trt_raw=tr
                ))
    return rows


# ---------------------------------------------------------------------------
# Differential metrics
# ---------------------------------------------------------------------------

def _get_fisher():
    try:
        from scipy.stats import fisher_exact
        return fisher_exact
    except ImportError:
        return None


def bh_adjust(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjustment with monotonicity enforcement."""
    n = len(pvals)
    if n == 0:
        return []
    order  = sorted(range(n), key=lambda i: pvals[i])
    q      = [1.0] * n
    cummin = 1.0
    for k in range(n - 1, -1, -1):
        i    = order[k]
        rank = k + 1
        val  = pvals[i] * n / rank
        cummin = min(cummin, val)
        q[i] = min(cummin, 1.0)
    return q


def finalize_diff_metrics(
    rows:           list[DiffRow],
    ctrl_total:     int,
    trt_total:      int,
    ctrl_total_mod: dict[str, int],
    trt_total_mod:  dict[str, int],
    pseudocount:    float,
    do_test:        bool,
) -> bool:
    """
    Fill cpm, log2FC and (optionally) Fisher p/q on each DiffRow.
    Returns True if a statistical test was actually run.
    """
    c_scale = 1e6 / ctrl_total if ctrl_total else 0.0
    t_scale = 1e6 / trt_total  if trt_total  else 0.0

    for r in rows:
        r.ctrl_cpm = r.ctrl_raw * c_scale
        r.trt_cpm  = r.trt_raw * t_scale
        r.log2fc   = math.log2((r.trt_cpm + pseudocount) / (r.ctrl_cpm + pseudocount))

    if not do_test:
        return False

    fisher = _get_fisher()
    if fisher is None:
        logger.warning(
            "scipy not available — skipping per-feature Fisher test. "
            "Install scipy or rerun with --no-test to silence this."
        )
        return False

    pvals: list[float] = []
    for r in rows:
        ct = ctrl_total_mod.get(r.mod_code, 0)
        tt = trt_total_mod.get(r.mod_code, 0)
        a  = r.ctrl_raw
        b  = max(ct - a, 0)
        c  = r.trt_raw
        d  = max(tt - c, 0)
        if (a + b) == 0 or (c + d) == 0:
            p = float("nan")
        else:
            try:
                _, p = fisher([[a, b], [c, d]])
            except Exception:
                p = float("nan")
        r.pvalue = p
        pvals.append(p)

    clean = [p if p == p else 1.0 for p in pvals]
    qs    = bh_adjust(clean)
    for r, p, q in zip(rows, pvals, qs):
        r.qvalue = q if p == p else float("nan")
    return True


# ---------------------------------------------------------------------------
# Replicon length helpers (for plotting)
# ---------------------------------------------------------------------------

def load_sequence_regions(gff_path: Path) -> dict[str, int]:
    """Parse '##sequence-region <seqid> <start> <end>' pragmas."""
    lengths: dict[str, int] = {}
    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.startswith("##sequence-region"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    lengths[parts[1]] = int(parts[3])
                except ValueError:
                    continue
    return lengths


def compute_replicon_lengths(
    features: list[GffFeature],
    sites:    list[MethylSite],
    declared: dict[str, int],
) -> dict[str, int]:
    observed_max: dict[str, int] = defaultdict(int)
    for f in features:
        observed_max[f.seqid] = max(observed_max[f.seqid], f.end)
    for s in sites:
        observed_max[s.chrom] = max(observed_max[s.chrom], s.start + 1)

    lengths: dict[str, int] = {}
    for seqid, obs in observed_max.items():
        lengths[seqid] = declared.get(seqid, obs) or obs
    return lengths


def _coord_unit(max_len: int) -> tuple[float, str]:
    if max_len >= 1_000_000:
        return 1_000_000.0, "Mb"
    if max_len >= 10_000:
        return 1_000.0, "kb"
    return 1.0, "bp"


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _output_stem_and_ext(output_name: str) -> tuple[str, str]:
    p   = Path(output_name)
    ext = p.suffix if p.suffix else ".tsv"
    return p.stem, ext


def _safe_label(label: str) -> str:
    return "".join(ch if (ch.isalnum() or ch in "-_.") else "_" for ch in label)


def write_condition_outputs(
    rows:       list[OutputRow],
    mod_codes:  set[str],
    output_dir: Path,
    stem:       str,
    ext:        str,
    label:      str,
    split:      bool,
) -> None:
    """Per-condition single listing (one set of file(s) per BED/GFF input)."""
    safe = _safe_label(label)
    if split:
        for mod in sorted(mod_codes):
            ml       = MOD_CODE_TO_LABEL.get(mod, mod)
            fpath    = output_dir / f"{stem}_{safe}_{ml}{ext}"
            mod_rows = [r for r in rows if r.mod_code == mod]
            with open(fpath, "w", encoding="utf-8") as fh:
                fh.write(OUTPUT_HEADER + "\n")
                for row in mod_rows:
                    fh.write(row.to_tsv() + "\n")
            logger.info("  [%s | %s]  %d row(s)  →  %s",
                        label, ml, len(mod_rows), fpath.name)
    else:
        fpath = output_dir / f"{stem}_{safe}{ext}"
        with open(fpath, "w", encoding="utf-8") as fh:
            fh.write(OUTPUT_HEADER + "\n")
            for row in rows:
                fh.write(row.to_tsv() + "\n")
        logger.info("  [%s | combined]  %d row(s)  →  %s",
                    label, len(rows), fpath.name)


def write_diff_outputs(
    rows:       list[DiffRow],
    mod_codes:  set[str],
    output_dir: Path,
    stem:       str,
    ext:        str,
    labels:     tuple[str, str],
    split:      bool,
    do_test:    bool,
) -> None:
    """Crossed differential table."""
    header = diff_header(labels, do_test)
    if split:
        for mod in sorted(mod_codes):
            ml       = MOD_CODE_TO_LABEL.get(mod, mod)
            fpath    = output_dir / f"{stem}_differential_{ml}{ext}"
            mod_rows = [r for r in rows if r.mod_code == mod]
            with open(fpath, "w", encoding="utf-8") as fh:
                fh.write(header + "\n")
                for row in mod_rows:
                    fh.write(row.to_tsv(do_test) + "\n")
            logger.info("  [differential | %s]  %d row(s)  →  %s",
                        ml, len(mod_rows), fpath.name)
    else:
        fpath = output_dir / f"{stem}_differential{ext}"
        with open(fpath, "w", encoding="utf-8") as fh:
            fh.write(header + "\n")
            for row in rows:
                fh.write(row.to_tsv(do_test) + "\n")
        logger.info("  [differential | combined]  %d row(s)  →  %s",
                    len(rows), fpath.name)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _lazy_import_plotting():
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        logger.error(
            "Plotting requires numpy and matplotlib (%s). "
            "Install them or rerun with --no-plot.", exc,
        )
        return None
    return np, plt


def _group_positions(sites: list[MethylSite]) -> dict[tuple[str, str], list[int]]:
    d: dict[tuple[str, str], list[int]] = defaultdict(list)
    for s in sites:
        d[(s.chrom, s.mod)].append(s.start)
    return d


def _plot_replicons(group_c, group_t, replicon_lengths) -> list[str]:
    reps = ({c for (c, _m) in group_c} | {c for (c, _m) in group_t})
    return sorted(reps & set(replicon_lengths))


def plot_along_replicons_overlay(
    ctrl_sites, trt_sites, replicon_lengths, mod_codes, totals, labels,
    output_dir, stem, fmt, bins, dpi,
):
    imported = _lazy_import_plotting()
    if imported is None:
        return None
    np, plt = imported

    ctrl_total, trt_total = totals
    gc, gt = _group_positions(ctrl_sites), _group_positions(trt_sites)
    replicons = _plot_replicons(gc, gt, replicon_lengths)
    if not replicons:
        logger.warning("Overlay plot: no replicon with sites — skipped.")
        return None

    mods = sorted(mod_codes)
    max_len = max(replicon_lengths[r] for r in replicons)
    scale, unit = _coord_unit(max_len)
    c_scale = 1e6 / ctrl_total if ctrl_total else 0.0
    t_scale = 1e6 / trt_total  if trt_total  else 0.0
    col_c, col_t = CONDITION_PALETTE

    nrows, ncols = len(replicons), len(mods)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(max(4.5 * ncols, 6), max(2.2 * nrows, 3)),
        squeeze=False, sharex="row",
    )

    for i, rep in enumerate(replicons):
        length  = replicon_lengths[rep]
        edges   = np.linspace(0, length, bins + 1)
        centres = (edges[:-1] + edges[1:]) / 2.0 / scale
        for j, mod in enumerate(mods):
            ax = axes[i][j]
            cc, _ = np.histogram(np.asarray(gc.get((rep, mod), []), dtype=float), bins=edges)
            tc, _ = np.histogram(np.asarray(gt.get((rep, mod), []), dtype=float), bins=edges)
            ax.plot(centres, cc * c_scale, drawstyle="steps-mid",
                    color=col_c, linewidth=0.9, label=labels[0])
            ax.plot(centres, tc * t_scale, drawstyle="steps-mid",
                    color=col_t, linewidth=0.9, label=labels[1])
            ax.fill_between(centres, cc * c_scale, step="mid", color=col_c, alpha=0.20)
            ax.fill_between(centres, tc * t_scale, step="mid", color=col_t, alpha=0.20)

            if i == 0:
                ax.set_title(f"{MOD_CODE_TO_LABEL.get(mod, mod)} ({mod})", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{rep}\nCPM / bin", fontsize=8)
            if i == nrows - 1:
                ax.set_xlabel(f"position ({unit})", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.margins(x=0)
            if i == 0 and j == ncols - 1:
                ax.legend(fontsize=8, frameon=True)

    fig.suptitle(
        "Methylation distribution along replicon(s): "
        f"{labels[0]} vs {labels[1]} (library-size normalised)",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fpath = output_dir / f"{stem}_along_replicons_overlay.{fmt}"
    fig.savefig(fpath, dpi=dpi)
    plt.close(fig)
    logger.info("  [plot] along-replicon overlay  →  %s", fpath.name)
    return fpath


def plot_along_replicons_diff(
    ctrl_sites, trt_sites, replicon_lengths, mod_codes, totals, labels,
    output_dir, stem, fmt, bins, dpi,
):
    imported = _lazy_import_plotting()
    if imported is None:
        return None
    np, plt = imported

    ctrl_total, trt_total = totals
    gc, gt = _group_positions(ctrl_sites), _group_positions(trt_sites)
    replicons = _plot_replicons(gc, gt, replicon_lengths)
    if not replicons:
        logger.warning("Difference plot: no replicon with sites — skipped.")
        return None

    mods = sorted(mod_codes)
    max_len = max(replicon_lengths[r] for r in replicons)
    scale, unit = _coord_unit(max_len)
    c_scale = 1e6 / ctrl_total if ctrl_total else 0.0
    t_scale = 1e6 / trt_total  if trt_total  else 0.0
    gain_col, loss_col = "#1a9850", "#762a83"   # gain in treatment / loss

    nrows, ncols = len(replicons), len(mods)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(max(4.5 * ncols, 6), max(2.2 * nrows, 3)),
        squeeze=False, sharex="row",
    )

    for i, rep in enumerate(replicons):
        length  = replicon_lengths[rep]
        edges   = np.linspace(0, length, bins + 1)
        centres = (edges[:-1] + edges[1:]) / 2.0 / scale
        for j, mod in enumerate(mods):
            ax = axes[i][j]
            cc, _ = np.histogram(np.asarray(gc.get((rep, mod), []), dtype=float), bins=edges)
            tc, _ = np.histogram(np.asarray(gt.get((rep, mod), []), dtype=float), bins=edges)
            delta = tc * t_scale - cc * c_scale
            ax.fill_between(centres, delta, 0, where=(delta >= 0), step="mid",
                            color=gain_col, alpha=0.85, linewidth=0)
            ax.fill_between(centres, delta, 0, where=(delta < 0), step="mid",
                            color=loss_col, alpha=0.85, linewidth=0)
            ax.axhline(0, color="black", linewidth=0.5)

            if i == 0:
                ax.set_title(f"{MOD_CODE_TO_LABEL.get(mod, mod)} ({mod})", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{rep}\nΔCPM / bin", fontsize=8)
            if i == nrows - 1:
                ax.set_xlabel(f"position ({unit})", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.margins(x=0)

    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=gain_col, alpha=0.85, label=f"gain in {labels[1]}"),
        Patch(facecolor=loss_col, alpha=0.85, label=f"loss in {labels[1]}"),
    ]
    fig.legend(handles=handles, loc="upper right", fontsize=8, frameon=True)
    fig.suptitle(
        f"Differential methylation along replicon(s): {labels[1]} − {labels[0]} (ΔCPM)",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fpath = output_dir / f"{stem}_along_replicons_diff.{fmt}"
    fig.savefig(fpath, dpi=dpi)
    plt.close(fig)
    logger.info("  [plot] along-replicon difference  →  %s", fpath.name)
    return fpath


def plot_per_replicon_feature(
    diff_rows, mod_codes, totals, labels, output_dir, stem, fmt, dpi,
):
    imported = _lazy_import_plotting()
    if imported is None:
        return None
    np, plt = imported

    ctrl_total, trt_total = totals
    c_scale = 1e6 / ctrl_total if ctrl_total else 0.0
    t_scale = 1e6 / trt_total  if trt_total  else 0.0

    # (replicon, feat_type, mod) -> [ctrl_raw, trt_raw]
    agg: dict[tuple[str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    replicons: set[str] = set()
    feat_types: set[str] = set()
    for r in diff_rows:
        key = (r.feature.seqid, r.feature.feat_type, r.mod_code)
        agg[key][0] += r.ctrl_raw
        agg[key][1] += r.trt_raw
        replicons.add(r.feature.seqid)
        feat_types.add(r.feature.feat_type)

    if not agg:
        logger.warning("Per-feature plot: no methylated features — skipped.")
        return None

    replicons_s = sorted(replicons)
    feats_s     = sorted(feat_types)
    mods        = sorted(mod_codes)
    col_c, col_t = CONDITION_PALETTE

    nrows, ncols = len(replicons_s), len(mods)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(max(1.6 * len(feats_s) * ncols, 7), 3.0 * nrows),
        squeeze=False, sharex="col",
    )

    x     = np.arange(len(feats_s))
    width = 0.4

    for i, rep in enumerate(replicons_s):
        for j, mod in enumerate(mods):
            ax = axes[i][j]
            ch = [agg.get((rep, ft, mod), [0, 0])[0] * c_scale for ft in feats_s]
            th = [agg.get((rep, ft, mod), [0, 0])[1] * t_scale for ft in feats_s]
            ax.bar(x - width / 2, ch, width=width, color=col_c, label=labels[0])
            ax.bar(x + width / 2, th, width=width, color=col_t, label=labels[1])

            if i == 0:
                ax.set_title(f"{MOD_CODE_TO_LABEL.get(mod, mod)} ({mod})", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{rep}\nCPM", fontsize=8)
            ax.set_xticks(x)
            ax.set_xticklabels(feats_s, rotation=30, ha="right", fontsize=8)
            ax.tick_params(axis="y", labelsize=7)
            ax.margins(x=0.02)
            if i == 0 and j == ncols - 1:
                ax.legend(fontsize=8, frameon=True)

    fig.suptitle(
        "Methylation per replicon and per feature: "
        f"{labels[0]} vs {labels[1]} (library-size normalised)",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fpath = output_dir / f"{stem}_per_replicon_feature.{fmt}"
    fig.savefig(fpath, dpi=dpi)
    plt.close(fig)
    logger.info("  [plot] per-replicon × feature comparison  →  %s", fpath.name)
    return fpath


def make_plots(
    features, ctrl_sites, trt_sites, diff_rows, mod_codes, totals, labels,
    declared_lengths, output_dir, output_name, fmt, bins, dpi,
) -> None:
    stem, _ext = _output_stem_and_ext(output_name)
    all_sites  = ctrl_sites + trt_sites
    replicon_lengths = compute_replicon_lengths(features, all_sites, declared_lengths)

    logger.info("Generating figure(s):")
    plot_along_replicons_overlay(
        ctrl_sites, trt_sites, replicon_lengths, mod_codes, totals, labels,
        output_dir, stem, fmt, bins, dpi,
    )
    plot_along_replicons_diff(
        ctrl_sites, trt_sites, replicon_lengths, mod_codes, totals, labels,
        output_dir, stem, fmt, bins, dpi,
    )

    distinct_feat_types = {r.feature.feat_type for r in diff_rows}
    if len(distinct_feat_types) > 1:
        plot_per_replicon_feature(
            diff_rows, mod_codes, totals, labels, output_dir, stem, fmt, dpi,
        )
    else:
        logger.info(
            "  [plot] per-feature plot skipped — only %d feature type with "
            "methylation.", len(distinct_feat_types),
        )


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(
    features, diff_rows, mod_codes, labels, totals, sources, did_test,
) -> None:
    sep, sep2 = "=" * 76, "-" * 76
    c_lab, t_lab = labels
    ctrl_total, trt_total = totals

    replicon_feat: dict[str, int] = defaultdict(int)
    feat_type_cnt: dict[str, int] = defaultdict(int)
    for f in features:
        replicon_feat[f.seqid]     += 1
        feat_type_cnt[f.feat_type] += 1

    # per (replicon, mod) summed counts for each condition
    rep_mod_c: dict[tuple[str, str], int] = defaultdict(int)
    rep_mod_t: dict[tuple[str, str], int] = defaultdict(int)
    mod_c: dict[str, int] = defaultdict(int)
    mod_t: dict[str, int] = defaultdict(int)
    for r in diff_rows:
        rep_mod_c[(r.feature.seqid, r.mod_code)] += r.ctrl_raw
        rep_mod_t[(r.feature.seqid, r.mod_code)] += r.trt_raw
        mod_c[r.mod_code] += r.ctrl_raw
        mod_t[r.mod_code] += r.trt_raw

    print(f"\n{sep}")
    print("  Differential Methylation × Annotation Cross-Reference — Summary")
    print(f"  {c_lab} source: {sources[0]}")
    print(f"  {t_lab} source: {sources[1]}")
    print(f"  Library sizes (total calls): {c_lab}={ctrl_total:,}  {t_lab}={trt_total:,}")
    print(sep)

    print(f"\n  Replicons found: {len(replicon_feat)}")
    print(f"  {sep2}")
    print(f"  {'Replicon':<36} {'Features':>10}")
    print(f"  {sep2}")
    for rep, cnt in sorted(replicon_feat.items()):
        print(f"  {rep:<36} {cnt:>10,}")

    print(f"\n  Feature types loaded: {len(feat_type_cnt)}")
    print(f"  {sep2}")
    print(f"  {'Type':<25} {'Count':>10}")
    print(f"  {sep2}")
    for ftype, cnt in sorted(feat_type_cnt.items(), key=lambda x: -x[1]):
        print(f"  {ftype:<25} {cnt:>10,}")

    print(f"\n  In-feature methylation by type (raw counts):")
    print(f"  {sep2}")
    print(f"  {'Code':<5} {'Label':<7} {c_lab:>14} {t_lab:>14} {'delta':>12}")
    print(f"  {sep2}")
    for mod in sorted(mod_codes):
        ml = MOD_CODE_TO_LABEL.get(mod, mod)
        cv, tv = mod_c.get(mod, 0), mod_t.get(mod, 0)
        print(f"  {mod:<5} {ml:<7} {cv:>14,} {tv:>14,} {tv - cv:>+12,}")

    print(f"\n  In-feature counts per replicon × mod ({c_lab} | {t_lab} | Δ):")
    print(f"  {sep2}")
    for rep in sorted(replicon_feat):
        print(f"  {rep}")
        for mod in sorted(mod_codes):
            ml = MOD_CODE_TO_LABEL.get(mod, mod)
            cv = rep_mod_c.get((rep, mod), 0)
            tv = rep_mod_t.get((rep, mod), 0)
            print(f"    {ml:<7} {cv:>12,} | {tv:>12,} | {tv - cv:>+12,}")

    if did_test:
        sig = [r for r in diff_rows if r.qvalue == r.qvalue and r.qvalue < 0.05]
        up   = [r for r in sig if r.log2fc > 0]
        down = [r for r in sig if r.log2fc < 0]
        print(f"\n  Differential features (Fisher exact, BH q < 0.05): {len(sig)}")
        print(f"    up in {t_lab}: {len(up)}   |   down in {t_lab}: {len(down)}")
        top = sorted(sig, key=lambda r: r.qvalue)[:10]
        if top:
            print(f"  {sep2}")
            print(f"  {'locus_tag':<14} {'gene':<8} {'mod':<5} "
                  f"{'log2FC':>8} {'qvalue':>10}")
            print(f"  {sep2}")
            for r in top:
                ml = MOD_CODE_TO_LABEL.get(r.mod_code, r.mod_code)
                print(f"  {r.feature.locus_tag:<14} {r.feature.gene:<8} {ml:<5} "
                      f"{r.log2fc:>8.3f} {r.qvalue:>10.2e}")
    else:
        print("\n  Per-feature statistical test: not run "
              "(scipy missing or --no-test).")

    print(f"\n{sep}\n")


# ---------------------------------------------------------------------------
# Seqid reconciliation
# ---------------------------------------------------------------------------

def build_seqid_map(features: list[GffFeature], sites: list[MethylSite]) -> dict[str, str]:
    gff_seqids = sorted({f.seqid for f in features})
    bed_seqids = sorted({s.chrom for s in sites})

    common = set(gff_seqids) & set(bed_seqids)
    if common:
        return {s: s for s in bed_seqids}

    if len(gff_seqids) != len(bed_seqids):
        logger.warning(
            "Chromosome name mismatch AND different counts "
            "(%d in GFF vs %d in methylation) — cannot auto-map. "
            "Use --seqid-map to supply an explicit TSV mapping.",
            len(gff_seqids), len(bed_seqids),
        )
        return {s: s for s in bed_seqids}

    mapping: dict[str, str] = {}
    for bed_name, gff_name in zip(bed_seqids, gff_seqids):
        mapping[bed_name] = gff_name
        logger.warning("Auto seqid mapping: '%s' → '%s'", bed_name, gff_name)
    return mapping


def apply_seqid_map(sites: list[MethylSite], seqid_map: dict[str, str]) -> list[MethylSite]:
    return [
        MethylSite(
            chrom  = seqid_map.get(s.chrom, s.chrom),
            start  = s.start,
            mod    = s.mod,
            strand = s.strand,
        )
        for s in sites
    ]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Differential cross-reference of a Prokka GFF3 annotation with "
            "methylation from two conditions (control vs treatment).\n\n"
            "Provide exactly one control input (-b/-m) and exactly one "
            "treatment input (-B/-M)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument("-g", "--gff", required=True, metavar="FILE",
                     help="Prokka GFF3 annotation file (shared by both conditions).")
    req.add_argument("-o", "--output-dir", required=True, metavar="DIR",
                     help="Directory where outputs and figures are written.")

    ctrl = parser.add_argument_group("control methylation input (exactly one)")
    ctrl_ex = ctrl.add_mutually_exclusive_group(required=True)
    ctrl_ex.add_argument("-b", "--bed", metavar="FILE",
                         help="Filtered 6-column BED (control).")
    ctrl_ex.add_argument("-m", "--methyl-gff", metavar="FILE",
                         help="Basecalling GFF3 (control).")

    trt = parser.add_argument_group("treatment methylation input (exactly one)")
    trt_ex = trt.add_mutually_exclusive_group(required=True)
    trt_ex.add_argument("-B", "--bed2", metavar="FILE",
                        help="Filtered 6-column BED (treatment).")
    trt_ex.add_argument("-M", "--methyl-gff2", metavar="FILE",
                        help="Basecalling GFF3 (treatment).")

    out = parser.add_argument_group("output options")
    out.add_argument("-n", "--name", default=None, metavar="FILENAME",
                     help="Base name for outputs. Default "
                          "'<gff_stem>_methylation.tsv'.")
    out.add_argument("--labels", nargs=2, default=["control", "treatment"],
                     metavar=("CONTROL", "TREATMENT"),
                     help="Condition labels used in columns/filenames/plots.")
    out.add_argument("--no-split", action="store_true",
                     help="Write single combined TSVs instead of one per mod type.")

    diff = parser.add_argument_group("differential options")
    diff.add_argument("--pseudocount", type=float, default=1.0, metavar="X",
                      help="Pseudocount (CPM units) added before log2FC. Default 1.0.")
    diff.add_argument("--no-test", action="store_true",
                      help="Disable the per-feature Fisher exact test + BH correction.")

    map_grp = parser.add_argument_group("seqid mapping")
    map_grp.add_argument("--seqid-map", default=None, metavar="FILE",
                         help="Two-column TSV mapping methylation seqids to GFF "
                              "seqids (applied to both conditions).")

    flt = parser.add_argument_group("filtering options")
    flt.add_argument("--features", nargs="+", default=None, metavar="TYPE",
                     help="GFF feature types to include (default: all non-meta).")
    flt.add_argument("--ignore-strand", action="store_true",
                     help="Count methylation on both strands per feature.")

    plot = parser.add_argument_group("plotting options")
    plot.add_argument("--no-plot", action="store_true",
                      help="Disable figure generation.")
    plot.add_argument("--plot-format", default="png",
                      choices=["png", "pdf", "svg"], metavar="FMT",
                      help="Figure format. Default png.")
    plot.add_argument("--bins", type=int, default=200, metavar="N",
                      help="Positional bins for along-replicon plots. Default 200.")
    plot.add_argument("--plot-dpi", type=int, default=150, metavar="DPI",
                      help="Raster resolution for png. Default 150.")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _load_explicit_map(path: str) -> dict[str, str]:
    m: dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m


def main() -> None:
    args = parse_args()

    gff_path   = Path(args.gff).resolve()
    output_dir = Path(args.output_dir).resolve()
    labels     = (str(args.labels[0]), str(args.labels[1]))

    if labels[0] == labels[1]:
        logger.error("The two condition labels must differ (got '%s' twice).",
                     labels[0])
        sys.exit(1)

    if not gff_path.is_file():
        logger.error("Annotation GFF file not found: %s", gff_path)
        sys.exit(1)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_name   = args.name or f"{gff_path.stem}_methylation.tsv"
    stem, ext     = _output_stem_and_ext(output_name)
    feature_types = set(args.features) if args.features else None
    split         = not args.no_split
    do_test       = not args.no_test

    ctrl_path, ctrl_src, ctrl_loader = loader_for(args.bed, args.methyl_gff)
    trt_path,  trt_src,  trt_loader  = loader_for(args.bed2, args.methyl_gff2)

    for p in (ctrl_path, trt_path):
        if not p.is_file():
            logger.error("Methylation input file not found: %s", p)
            sys.exit(1)

    explicit_map = _load_explicit_map(args.seqid_map) if args.seqid_map else None
    if explicit_map is not None:
        logger.info("Loaded %d explicit seqid mapping(s).", len(explicit_map))

    logger.info("Annotation GFF     : %s", gff_path)
    logger.info("Control (%s)  : %s", labels[0], ctrl_src)
    logger.info("Treatment (%s): %s", labels[1], trt_src)
    logger.info("Output directory   : %s", output_dir)
    logger.info("Output base name   : %s", output_name)
    logger.info("Split by mod type  : %s", split)
    logger.info("Per-feature test   : %s", "Fisher+BH" if do_test else "disabled")
    logger.info(
        "Feature types      : %s",
        ", ".join(sorted(feature_types)) if feature_types
        else f"all (excluding {', '.join(sorted(GFF_META_TYPES))})",
    )
    logger.info("Strand-aware       : %s", not args.ignore_strand)

    features = load_annotation_gff(gff_path, feature_types=feature_types)
    if not features:
        logger.error("No features loaded from annotation GFF.")
        sys.exit(1)

    ctrl_sites = ctrl_loader(ctrl_path)
    trt_sites  = trt_loader(trt_path)
    if not ctrl_sites or not trt_sites:
        logger.error("A condition has no methylation sites — cannot run differential.")
        sys.exit(1)

    # Seqid reconciliation (shared map from features + both conditions)
    if explicit_map is not None:
        seqid_map = explicit_map
    else:
        seqid_map = build_seqid_map(features, ctrl_sites + trt_sites)
    ctrl_sites = apply_seqid_map(ctrl_sites, seqid_map)
    trt_sites  = apply_seqid_map(trt_sites, seqid_map)

    gff_seqids = {f.seqid for f in features}
    meth_seqids = {s.chrom for s in ctrl_sites} | {s.chrom for s in trt_sites}
    overlap = gff_seqids & meth_seqids
    if not overlap:
        logger.error(
            "No chromosome names overlap between annotation (%s …) and "
            "methylation (%s …) after mapping. Supply --seqid-map.",
            ", ".join(sorted(gff_seqids)[:3]),
            ", ".join(sorted(meth_seqids)[:3]),
        )
        sys.exit(1)
    logger.info("Seqid overlap after mapping: %d replicon(s) in common.", len(overlap))

    mod_codes  = {s.mod for s in ctrl_sites} | {s.mod for s in trt_sites}
    idx_ctrl   = build_methyl_index(ctrl_sites)
    idx_trt    = build_methyl_index(trt_sites)

    ctrl_total = len(ctrl_sites)
    trt_total  = len(trt_sites)
    ctrl_total_mod = dict(Counter(s.mod for s in ctrl_sites))
    trt_total_mod  = dict(Counter(s.mod for s in trt_sites))
    totals = (ctrl_total, trt_total)

    # Per-condition single listings
    logger.info("Writing per-condition listing(s):")
    ctrl_rows = cross_reference(features, idx_ctrl, mod_codes, args.ignore_strand)
    trt_rows  = cross_reference(features, idx_trt,  mod_codes, args.ignore_strand)
    write_condition_outputs(ctrl_rows, mod_codes, output_dir, stem, ext, labels[0], split)
    write_condition_outputs(trt_rows,  mod_codes, output_dir, stem, ext, labels[1], split)

    # Crossed differential table
    logger.info("Building crossed differential table…")
    diff_rows = cross_reference_diff(features, idx_ctrl, idx_trt, mod_codes,
                                     args.ignore_strand)
    did_test = finalize_diff_metrics(
        diff_rows, ctrl_total, trt_total, ctrl_total_mod, trt_total_mod,
        args.pseudocount, do_test,
    )
    logger.info("Writing differential output(s):")
    write_diff_outputs(diff_rows, mod_codes, output_dir, stem, ext, labels,
                       split, did_test)

    # Figures
    if not args.no_plot:
        declared_lengths = load_sequence_regions(gff_path)
        make_plots(
            features, ctrl_sites, trt_sites, diff_rows, mod_codes, totals, labels,
            declared_lengths, output_dir, output_name,
            args.plot_format, args.bins, args.plot_dpi,
        )
    else:
        logger.info("Plotting disabled (--no-plot).")

    print_summary(features, diff_rows, mod_codes, labels, totals,
                  (ctrl_src, trt_src), did_test)


if __name__ == "__main__":
    main()
