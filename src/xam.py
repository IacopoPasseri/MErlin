#!/usr/bin/env python3
"""
xam.py
---------------------------------------------------------------------------------
Cross-references a GFF3 genome annotation with methylation basecalling.

Methylation input:
    -b / --bed          Filtered 6-column BED file from mbr.py
    -m / --methyl-gff   GFF3 file produced by the basecalling where every
                        feature record represents a single methylation call.

For every annotated feature the script reports:
    - replicon (chromosome / contig)
    - feature type  (CDS, tRNA, rRNA, …)
    - strand
    - feature length (bp)
    - methylation type (e.g. m = 5mC, a = 6mA …)
    - raw count    : number of methylated base positions inside the feature
                     on the same strand
    - norm count   : raw_count / feature_length_bp
                     (fraction of bases in the feature that are methylated)

---------------------------------------------------------------------------------
Overlap rule
---------------------------------------------------------------------------------
A methylation site at position P (0-based) is counted for a feature
[feat_start, feat_end) if and only if:
    1.  feat_start <= P < feat_end          (position is inside the feature)
    2.  site strand == feature strand       (same strand; override: --ignore-strand)
    3.  site chrom  == feature seqid        (same replicon — always enforced)

---------------------------------------------------------------------------------
Basecalling GFF format
---------------------------------------------------------------------------------
Each data line in the methylation GFF is expected to follow GFF3 conventions:
    col 0  seqid      — replicon name
    col 2  type       — modification type identifier; the script looks for:
                          • a known mod code directly  ("m", "a", "h", …)
                          • a known mod label directly  ("5mC", "6mA", …)
                          • common synonyms: "modified_base", "methylation",
                            "CpG", "dam", "dcm", "GATC", "CCWGG"
                          • fallback: the raw type string is kept as-is so
                            that novel types are preserved rather than dropped.
                        The mod code is also looked up in the 'Name', 'mod',
                        'modification', and 'type' GFF attributes as a
                        secondary source.
    col 3  start      — 1-based inclusive start (converted to 0-based)
    col 6  strand     — '+' or '-'

---------------------------------------------------------------------------------
Usage
---------------------------------------------------------------------------------
    # BED input (original behaviour)
    python3 xam.py \\
        -g annotation.gff -b filtered.bed -o results/

    # GFF methylation input
    python3 xam.py \\
        -g annotation.gff -m basecalling.gff -o results/

    # Custom output name
    python3 xam.py \\
        -g annotation.gff -m basecalling.gff -o results/ -n my_results.tsv

    # Single combined file
    python3 xam.py \\
        -g annotation.gff -b filtered.bed -o results/ --no-split

    # Only CDS and tRNA features
    python3 xam.py \\
        -g annotation.gff -b filtered.bed -o results/ --features CDS tRNA

    # Count both strands per feature
    python3 xam.py \\
        -g annotation.gff -b filtered.bed -o results/ --ignore-strand

"""

import argparse
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
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

# Modification code → short label (used in output filenames)
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

# Short label → human-readable description
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

# Reverse map: short label → mod code
MOD_LABEL_TO_CODE: dict[str, str] = {v: k for k, v in MOD_CODE_TO_LABEL.items()}

# Common synonyms found in basecalling GFF type fields → mod code
METHYL_GFF_TYPE_SYNONYMS: dict[str, str] = {
    # CpG / cytosine methylation
    "cpg":           "m",
    "5mc":           "m",
    "m5C":           "m",
    "5-mc":          "m",
    "methylation":   "m",   # generic fallback to 5mC
    "modified_base": "m",
    # Dam / adenine methylation
    "6ma":           "a",
    "m6A":           "a",
    "6-ma":          "a",
    "dam":           "a",
    "gatc":          "a",
    # Dcm / cytosine methylation (CCWGG context)
    "dcm":           "m",
    "ccwgg":         "m",
    # hydroxymethyl
    "5hmc":          "h",
    "5-hmc":         "h",
}

# GFF3 column indices (0-based)
GFF_SEQID  = 0
GFF_TYPE   = 2
GFF_START  = 3   # 1-based inclusive
GFF_END    = 4   # 1-based inclusive
GFF_STRAND = 6
GFF_ATTRS  = 8

# BED column indices (0-based)
BED_CHROM  = 0
BED_START  = 1   # 0-based
BED_MOD    = 3
BED_STRAND = 5

# GFF meta-feature types to exclude from annotation loading
GFF_META_TYPES = {"region", "sequence_feature", "sequence_alteration"}

OUTPUT_HEADER = "\t".join([
    "replicon",
    "feature_type",
    "feat_start",       # 0-based
    "feat_end",         # 0-based exclusive
    "strand",
    "feat_length_bp",
    "locus_tag",
    "gene",
    "product",
    "mod_code",
    "mod_label",
    "mod_description",
    "raw_count",
    "norm_count",       # raw_count / feat_length_bp
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
    strand:     str    # '+', '-', or '.'
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
    start:  int    # 0-based position of the modified base
    mod:    str    # normalised modification code (e.g. 'm', 'a')
    strand: str


@dataclass
class OutputRow:
    """One result row: one annotated feature × one mod type."""
    feature:   GffFeature
    mod_code:  str
    raw_count: int

    @property
    def norm_count(self) -> float:
        """raw_count / feature_length_bp (fraction of methylated bases)."""
        if self.feature.length == 0:
            return 0.0
        return self.raw_count / self.feature.length

    def to_tsv(self) -> str:
        f     = self.feature
        label = MOD_CODE_TO_LABEL.get(self.mod_code, self.mod_code)
        desc  = MOD_LABEL_TO_DESC.get(label, "unknown modification")
        return "\t".join([
            f.seqid,
            f.feat_type,
            str(f.start),
            str(f.end),
            f.strand,
            str(f.length),
            f.locus_tag,
            f.gene,
            f.product,
            self.mod_code,
            label,
            desc,
            str(self.raw_count),
            f"{self.norm_count:.8f}",
        ])


# ---------------------------------------------------------------------------
# Shared GFF3 attribute parser
# ---------------------------------------------------------------------------

def _parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse a GFF3 attribute string into a key → value dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            key, _, val = part.partition("=")
            attrs[key.strip()] = val.strip()
    return attrs


def _split_gff_line(raw: str) -> list[str] | None:
    """
    Split a GFF3 line into columns.
    Tries tab-split first (standard); falls back to whitespace if < 9 cols.
    Returns None if the line cannot be parsed into >= 9 columns.
    """
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
    """
    Load annotated features from a Prokka GFF3 file.

    Prokka embeds the genome FASTA at the end of the GFF3 after a '##FASTA'
    directive — parsing stops there automatically.

    Parameters
    ----------
    gff_path : Path
        Path to the Prokka .gff file.
    feature_types : set[str] | None
        If given, only features of these types are loaded.
        If None, all types except GFF_META_TYPES are loaded.

    Returns
    -------
    list[GffFeature]
    """
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
                start_0 = int(cols[GFF_START]) - 1   # GFF 1-based → 0-based
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
# BED methylation loader
# ---------------------------------------------------------------------------

def load_bed(bed_path: Path) -> list[MethylSite]:
    """
    Load methylation sites from the 6-column filtered BED file produced
    by filter_bedmethyl.py.

    Columns used:
        0  chrom
        1  start  (0-based)
        3  mod code
        5  strand

    Returns
    -------
    list[MethylSite]
    """
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


# ---------------------------------------------------------------------------
# Basecalling GFF methylation loader
# ---------------------------------------------------------------------------

def _resolve_mod_code(feat_type: str, attrs: dict[str, str]) -> str:
    """
    Determine the normalised modification code from a basecalling GFF feature.

    Resolution order:
        1. feat_type is already a known mod code  ("m", "a", …)
        2. feat_type is a known short label        ("5mC", "6mA", …)
        3. feat_type (lower-cased) is a known synonym
        4. Attributes 'mod', 'modification', 'Name', or 'type' contain a
           known code or label
        5. Return feat_type as-is (preserves novel/unknown types)

    Parameters
    ----------
    feat_type : str
        Value of column 2 in the basecalling GFF line.
    attrs : dict[str, str]
        Parsed GFF3 attributes from column 8.

    Returns
    -------
    str  — normalised mod code or the raw feat_type string if unresolved.
    """
    # 1. Direct code match
    if feat_type in MOD_CODE_TO_LABEL:
        return feat_type

    # 2. Short label match
    if feat_type in MOD_LABEL_TO_CODE:
        return MOD_LABEL_TO_CODE[feat_type]

    # 3. Synonym match (case-insensitive)
    lower = feat_type.lower()
    if lower in METHYL_GFF_TYPE_SYNONYMS:
        return METHYL_GFF_TYPE_SYNONYMS[lower]

    # 4. Attribute look-up
    for attr_key in ("mod", "modification", "Name", "type"):
        val = attrs.get(attr_key, "")
        if val in MOD_CODE_TO_LABEL:
            return val
        if val in MOD_LABEL_TO_CODE:
            return MOD_LABEL_TO_CODE[val]
        val_lower = val.lower()
        if val_lower in METHYL_GFF_TYPE_SYNONYMS:
            return METHYL_GFF_TYPE_SYNONYMS[val_lower]

    # 5. Fallback — return type string unchanged
    return feat_type


def load_methyl_gff(gff_path: Path) -> list[MethylSite]:
    """
    Load methylation sites from a basecalling GFF3 file.

    Each data line is interpreted as one methylation call:
        col 0  seqid   — replicon / chromosome
        col 2  type    — used to derive the mod code (see _resolve_mod_code)
        col 3  start   — 1-based inclusive (converted to 0-based)
        col 6  strand  — '+' or '-'

    Lines starting with '#' and the embedded '##FASTA' section are skipped.
    Lines where the type resolves to a meta-feature name (region, etc.) are
    also skipped, as they are almost certainly annotation lines mixed in by
    mistake.

    Parameters
    ----------
    gff_path : Path
        Path to the basecalling GFF3 file.

    Returns
    -------
    list[MethylSite]
    """
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

            # Skip obvious annotation meta-features
            if feat_type in GFF_META_TYPES:
                skipped_meta += 1
                continue

            attrs    = _parse_attributes(cols[GFF_ATTRS])
            mod_code = _resolve_mod_code(feat_type, attrs)

            # Warn once per unrecognised type (still kept in output)
            if (
                mod_code not in MOD_CODE_TO_LABEL
                and mod_code not in unknown_types
            ):
                unknown_types.add(mod_code)
                logger.warning(
                    "Methylation GFF: unrecognised type '%s' — kept as-is. "
                    "Add it to MOD_CODE_TO_LABEL / METHYL_GFF_TYPE_SYNONYMS "
                    "for a human-readable label.",
                    mod_code,
                )

            try:
                # GFF is 1-based inclusive; we need 0-based position
                start_0 = int(cols[GFF_START]) - 1
            except ValueError:
                logger.warning(
                    "Methylation GFF line %d: cannot parse start — skipped.", line_no
                )
                skipped_bad += 1
                continue

            strand = cols[GFF_STRAND]

            sites.append(MethylSite(
                chrom  = cols[GFF_SEQID],
                start  = start_0,
                mod    = mod_code,
                strand = strand,
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


# ---------------------------------------------------------------------------
# Methylation index
# ---------------------------------------------------------------------------

def build_methyl_index(
    sites: list[MethylSite],
) -> dict[tuple[str, str, str], list[int]]:
    """
    Build a lookup dict: (chrom, strand, mod_code) → sorted list of 0-based
    positions.  Enables O(log n) range queries per feature.
    """
    index: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    for site in sites:
        index[(site.chrom, site.strand, site.mod)].append(site.start)
    for key in index:
        index[key].sort()
    return index


# ---------------------------------------------------------------------------
# Overlap counting — strict position-based, binary-search accelerated
# ---------------------------------------------------------------------------

def count_overlaps(
    feature:       GffFeature,
    index:         dict[tuple[str, str, str], list[int]],
    all_mod_codes: set[str],
    ignore_strand: bool,
) -> dict[str, int]:
    """
    Count methylation positions inside a feature per mod type.

    A site at 0-based position P overlaps feature [start, end) when:
        feature.start <= P < feature.end

    Strand is matched unless ignore_strand=True or feature.strand == '.'.

    Returns
    -------
    dict: mod_code → count  (only codes with count > 0 included)
    """
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

            # Binary search for first position >= feature.start
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


# ---------------------------------------------------------------------------
# Cross-reference
# ---------------------------------------------------------------------------

def cross_reference(
    features:      list[GffFeature],
    index:         dict[tuple[str, str, str], list[int]],
    all_mod_codes: set[str],
    ignore_strand: bool,
) -> list[OutputRow]:
    """
    Return one OutputRow per (feature × mod_code) where raw_count > 0.
    Features with no overlapping methylation are excluded from the output.
    """
    rows: list[OutputRow] = []
    for feat in features:
        counts = count_overlaps(feat, index, all_mod_codes, ignore_strand)
        for mod in sorted(all_mod_codes):
            raw = counts.get(mod, 0)
            if raw > 0:
                rows.append(OutputRow(feature=feat, mod_code=mod, raw_count=raw))
    return rows


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _output_stem_and_ext(output_name: str) -> tuple[str, str]:
    p   = Path(output_name)
    ext = p.suffix if p.suffix else ".tsv"
    return p.stem, ext


def write_outputs(
    rows:        list[OutputRow],
    mod_codes:   set[str],
    output_dir:  Path,
    output_name: str,
    split:       bool,
) -> None:
    """
    Write result TSV(s).

    split=True  (default) → one file per mod type: <stem>_<MOD_LABEL><ext>
    split=False           → single combined file:   <output_name>
    """
    stem, ext = _output_stem_and_ext(output_name)

    if split:
        for mod in sorted(mod_codes):
            label    = MOD_CODE_TO_LABEL.get(mod, mod)
            fpath    = output_dir / f"{stem}_{label}{ext}"
            mod_rows = [r for r in rows if r.mod_code == mod]
            with open(fpath, "w", encoding="utf-8") as fh:
                fh.write(OUTPUT_HEADER + "\n")
                for row in mod_rows:
                    fh.write(row.to_tsv() + "\n")
            logger.info("  [%s]  %d row(s)  →  %s", label, len(mod_rows), fpath.name)
    else:
        fpath = output_dir / output_name
        with open(fpath, "w", encoding="utf-8") as fh:
            fh.write(OUTPUT_HEADER + "\n")
            for row in rows:
                fh.write(row.to_tsv() + "\n")
        logger.info("Combined output: %s  (%d rows)", fpath, len(rows))


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(
    features:      list[GffFeature],
    rows:          list[OutputRow],
    mod_codes:     set[str],
    methyl_source: str,
) -> None:
    sep  = "=" * 68
    sep2 = "-" * 68

    replicon_feat: dict[str, int]            = defaultdict(int)
    feat_type_cnt: dict[str, int]            = defaultdict(int)
    rep_mod_hits:  dict[tuple[str,str], int] = defaultdict(int)
    mod_total:     dict[str, int]            = defaultdict(int)

    for f in features:
        replicon_feat[f.seqid]     += 1
        feat_type_cnt[f.feat_type] += 1

    for row in rows:
        rep_mod_hits[(row.feature.seqid, row.mod_code)] += row.raw_count
        mod_total[row.mod_code]                         += row.raw_count

    print(f"\n{sep}")
    print("  Methylation × Annotation Cross-Reference — Summary")
    print(f"  Methylation source: {methyl_source}")
    print(sep)

    # Replicons
    print(f"\n  Replicons found: {len(replicon_feat)}")
    print(f"  {sep2}")
    print(f"  {'Replicon':<36} {'Features':>10}")
    print(f"  {sep2}")
    for rep, cnt in sorted(replicon_feat.items()):
        print(f"  {rep:<36} {cnt:>10,}")

    # Feature types
    print(f"\n  Feature types loaded: {len(feat_type_cnt)}")
    print(f"  {sep2}")
    print(f"  {'Type':<25} {'Count':>10}")
    print(f"  {sep2}")
    for ftype, cnt in sorted(feat_type_cnt.items(), key=lambda x: -x[1]):
        print(f"  {ftype:<25} {cnt:>10,}")

    # Methylation types
    print(f"\n  Methylation types: {len(mod_codes)}")
    print(f"  {sep2}")
    print(f"  {'Code':<6}  {'Label':<8}  {'Description':<30} {'Total hits':>12}")
    print(f"  {sep2}")
    for mod, total in sorted(mod_total.items(), key=lambda x: -x[1]):
        label = MOD_CODE_TO_LABEL.get(mod, mod)
        desc  = MOD_LABEL_TO_DESC.get(label, "unknown modification")
        print(f"  {mod:<6}  {label:<8}  {desc:<30} {total:>12,}")

    # Per-replicon × mod breakdown
    if mod_codes:
        sorted_mods = sorted(mod_codes)
        labels      = [MOD_CODE_TO_LABEL.get(m, m) for m in sorted_mods]
        header_mods = "  ".join(f"{lb:>10}" for lb in labels)
        print(f"\n  Methylation hits per replicon × mod type:")
        print(f"  {sep2}")
        print(f"  {'Replicon':<36}  {header_mods}")
        print(f"  {sep2}")
        for rep in sorted(replicon_feat):
            row_vals = "  ".join(
                f"{rep_mod_hits.get((rep, m), 0):>10,}" for m in sorted_mods
            )
            print(f"  {rep:<36}  {row_vals}")

    print(f"\n{sep}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Cross-reference a Prokka GFF3 annotation with methylation data\n"
            "to count methylation events per genomic feature.\n\n"
            "Methylation input: provide exactly one of -b/--bed or -m/--methyl-gff.\n\n"
            "By default one TSV is written per methylation type:\n"
            "  <stem>_5mC.tsv, <stem>_6mA.tsv, …"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  BED input (original):\n"
            "    python3 xam.py \\\n"
            "        -g annotation.gff -b filtered.bed -o results/\n\n"
            "  GFF methylation input:\n"
            "    python3 xam.py \\\n"
            "        -g annotation.gff -m basecalling.gff -o results/\n\n"
            "  Custom base name:\n"
            "    python3 xam.py \\\n"
            "        -g annotation.gff -m basecalling.gff -o results/ -n my_results.tsv\n"
            "    # produces: results/my_results_5mC.tsv, results/my_results_6mA.tsv …\n\n"
            "  Single combined file:\n"
            "    python3 xam.py \\\n"
            "        -g annotation.gff -b filtered.bed -o results/ --no-split\n\n"
            "  CDS and tRNA only:\n"
            "    python3 xam.py \\\n"
            "        -g annotation.gff -b filtered.bed -o results/ --features CDS tRNA\n\n"
            "  Ignore strand:\n"
            "    python3 xam.py \\\n"
            "        -g annotation.gff -b filtered.bed -o results/ --ignore-strand\n"
        ),
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-g", "--gff",
        required=True, metavar="FILE",
        help="Prokka GFF3 annotation file.",
    )
    req.add_argument(
        "-o", "--output-dir",
        required=True, metavar="DIR",
        help="Directory where output TSV file(s) will be saved.",
    )

    meth = parser.add_argument_group(
        "methylation input (provide exactly one)"
    )
    meth_ex = meth.add_mutually_exclusive_group(required=True)
    meth_ex.add_argument(
        "-b", "--bed",
        metavar="FILE",
        help=(
            "Filtered 6-column BED file from filter_bedmethyl.py "
            "(chrom, start, end, mod_code, score, strand)."
        ),
    )
    meth_ex.add_argument(
        "-m", "--methyl-gff",
        metavar="FILE",
        help=(
            "Basecalling GFF3 file where each feature represents one "
            "methylation call. The modification type is derived from "
            "column 2 (feature type) or GFF attributes."
        ),
    )

    out = parser.add_argument_group("output options")
    out.add_argument(
        "-n", "--name",
        default=None, metavar="FILENAME",
        help=(
            "Base name for the output file(s) "
            "(e.g. 'results.tsv' → 'results_5mC.tsv', 'results_6mA.tsv', …). "
            "Defaults to '<gff_stem>_methylation_by_feature.tsv'."
        ),
    )
    out.add_argument(
        "--no-split",
        action="store_true",
        help="Write a single combined TSV instead of one file per methylation type.",
    )

    map_grp = parser.add_argument_group("seqid mapping (use when BED and GFF use different chromosome names)")
    map_grp.add_argument(
        "--seqid-map",
        default=None, metavar="FILE",
        help=(
            "Two-column TSV file (no header) mapping BED chromosome names to GFF "
            "seqid names: <bed_name>\\t<gff_name>. "
            "If omitted, the script attempts auto-detection by matching "
            "lexicographically sorted names in order."
        ),
    )

    flt = parser.add_argument_group("filtering options")
    flt.add_argument(
        "--features",
        nargs="+", default=None, metavar="TYPE",
        help=(
            "GFF feature types to include (e.g. CDS tRNA rRNA gene). "
            "Default: all types except meta-features "
            f"({', '.join(sorted(GFF_META_TYPES))})."
        ),
    )
    flt.add_argument(
        "--ignore-strand",
        action="store_true",
        help=(
            "Count methylation on both strands for every feature. "
            "Default: only sites on the same strand as the feature are counted."
        ),
    )

    plot = parser.add_argument_group("plotting options")
    plot.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable figure generation (only TSV output is written).",
    )
    plot.add_argument(
        "--plot-format",
        default="png", choices=["png", "pdf", "svg"], metavar="FMT",
        help="Output format for the figures (png, pdf, svg). Default: png.",
    )
    plot.add_argument(
        "--bins",
        type=int, default=200, metavar="N",
        help=(
            "Number of positional bins used for the along-replicon "
            "distribution plot. Default: 200."
        ),
    )
    plot.add_argument(
        "--plot-dpi",
        type=int, default=150, metavar="DPI",
        help="Raster resolution for png output. Default: 150.",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def build_seqid_map(features: list[GffFeature], sites: list[MethylSite]) -> dict[str, str]:
    """
    Auto-detect a chromosome/contig name mapping from methylation seqids
    to annotation seqids when they differ (e.g. NCBI accessions in the BED
    vs. local contig names in the GFF3).

    Strategy: match by rank order — the Nth distinct name in the BED is
    assumed to correspond to the Nth distinct name in the GFF (both sorted
    lexicographically, which preserves numeric order for typical assemblies).

    A warning is emitted for every mapping applied so the user can verify.

    Returns
    -------
    dict : bed_seqid → gff_seqid   (identity entries included for unchanged names)
    """
    gff_seqids  = sorted({f.seqid for f in features})
    bed_seqids  = sorted({s.chrom for s in sites})

    # If there is already a full overlap, no mapping is needed
    common = set(gff_seqids) & set(bed_seqids)
    if common:
        return {s: s for s in bed_seqids}

    if len(gff_seqids) != len(bed_seqids):
        logger.warning(
            "Chromosome name mismatch AND different counts "
            "(%d in GFF vs %d in BED) — cannot auto-map. "
            "Use --seqid-map to supply an explicit TSV mapping.",
            len(gff_seqids), len(bed_seqids),
        )
        return {s: s for s in bed_seqids}

    mapping: dict[str, str] = {}
    for bed_name, gff_name in zip(bed_seqids, gff_seqids):
        mapping[bed_name] = gff_name
        logger.warning(
            "Auto seqid mapping: '%s' (BED) → '%s' (GFF)", bed_name, gff_name
        )
    return mapping


def apply_seqid_map(sites: list[MethylSite], seqid_map: dict[str, str]) -> list[MethylSite]:
    """Return a new list of MethylSite with chrom names remapped."""
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
# Plotting
# ---------------------------------------------------------------------------

def _lazy_import_plotting():
    """
    Import matplotlib (Agg backend) and numpy only when plotting is requested,
    so the cross-referencing core has no hard dependency on them.
    """
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment dependent
        logger.error(
            "Plotting requires numpy and matplotlib (%s). "
            "Install them or rerun with --no-plot.", exc,
        )
        return None
    return np, plt


def load_sequence_regions(gff_path: Path) -> dict[str, int]:
    """
    Parse '##sequence-region <seqid> <start> <end>' pragmas from a GFF3 file
    to recover the full length of every replicon.

    Returns
    -------
    dict : seqid → length (bp). Empty if no pragmas are present.
    """
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
    features:    list[GffFeature],
    sites:       list[MethylSite],
    declared:    dict[str, int],
) -> dict[str, int]:
    """
    Resolve a length (bp) for every replicon seen in features or sites.

    Priority: declared '##sequence-region' length; otherwise the maximum
    observed coordinate (feature end or methylation position + 1) is used as
    a conservative fallback so binning still spans the data.
    """
    observed_max: dict[str, int] = defaultdict(int)
    for f in features:
        observed_max[f.seqid] = max(observed_max[f.seqid], f.end)
    for s in sites:
        observed_max[s.chrom] = max(observed_max[s.chrom], s.start + 1)

    lengths: dict[str, int] = {}
    for seqid, obs in observed_max.items():
        lengths[seqid] = declared.get(seqid, obs) or obs
    return lengths


def _mod_color_map(mod_codes, plt):
    """Stable mod_code → colour mapping using the default tab10 cycle."""
    cmap   = plt.get_cmap("tab10")
    ordered = sorted(mod_codes)
    return {mod: cmap(i % 10) for i, mod in enumerate(ordered)}


def _coord_unit(max_len: int) -> tuple[float, str]:
    """Pick a sensible coordinate unit (bp / kb / Mb) for axis labelling."""
    if max_len >= 1_000_000:
        return 1_000_000.0, "Mb"
    if max_len >= 10_000:
        return 1_000.0, "kb"
    return 1.0, "bp"


def plot_along_replicons(
    sites:             list[MethylSite],
    replicon_lengths:  dict[str, int],
    mod_codes:         set[str],
    output_dir:        Path,
    stem:              str,
    fmt:               str,
    bins:              int,
    dpi:               int,
):
    """
    Methylation distribution ALONG each replicon, split by methylation type.

    Layout: a grid of small multiples with one row per replicon and one
    column per methylation type. Each cell is a positional histogram of
    methylation calls (both strands combined) across `bins` equal-width
    windows spanning the replicon.

    Returns the written figure path, or None if there is nothing to plot.
    """
    imported = _lazy_import_plotting()
    if imported is None:
        return None
    np, plt = imported

    # Group methylation positions by (replicon, mod_code)
    pos_by: dict[tuple[str, str], list[int]] = defaultdict(list)
    for s in sites:
        pos_by[(s.chrom, s.mod)].append(s.start)

    replicons = sorted(
        {chrom for (chrom, _mod) in pos_by}
        & set(replicon_lengths)
    )
    if not replicons:
        logger.warning("Along-replicon plot: no replicon with sites — skipped.")
        return None

    mods    = sorted(mod_codes)
    colours = _mod_color_map(mod_codes, plt)
    max_len = max(replicon_lengths[r] for r in replicons)
    scale, unit = _coord_unit(max_len)

    nrows, ncols = len(replicons), len(mods)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(max(4.5 * ncols, 6), max(2.2 * nrows, 3)),
        squeeze=False,
        sharex="row",
    )

    for i, rep in enumerate(replicons):
        length = replicon_lengths[rep]
        edges  = np.linspace(0, length, bins + 1)
        centres = (edges[:-1] + edges[1:]) / 2.0 / scale
        for j, mod in enumerate(mods):
            ax = axes[i][j]
            positions = pos_by.get((rep, mod), [])
            if positions:
                counts, _ = np.histogram(np.asarray(positions), bins=edges)
                ax.fill_between(
                    centres, counts, step="mid",
                    color=colours[mod], alpha=0.85, linewidth=0,
                )
                ax.plot(centres, counts, drawstyle="steps-mid",
                        color=colours[mod], linewidth=0.6)
            else:
                ax.text(0.5, 0.5, "no calls", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="grey")

            label = MOD_CODE_TO_LABEL.get(mod, mod)
            if i == 0:
                ax.set_title(f"{label} ({mod})", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{rep}\ncalls / bin", fontsize=8)
            if i == nrows - 1:
                ax.set_xlabel(f"position ({unit})", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.margins(x=0)

    fig.suptitle(
        "Methylation distribution along replicon(s), split by methylation type",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    fpath = output_dir / f"{stem}_along_replicons.{fmt}"
    fig.savefig(fpath, dpi=dpi)
    plt.close(fig)
    logger.info("  [plot] along-replicon distribution  →  %s", fpath.name)
    return fpath


def plot_per_replicon_per_feature(
    rows:        list[OutputRow],
    mod_codes:   set[str],
    output_dir:  Path,
    stem:        str,
    fmt:         str,
    dpi:         int,
):
    """
    Methylation distribution per replicon AND per feature type, split by
    methylation type.

    Layout: one subplot (row) per replicon; within each, grouped bars give
    the total number of methylated positions (summed raw_count) for every
    feature type, one bar group per feature type and one coloured bar per
    methylation type.

    Returns the written figure path, or None if there is nothing to plot.
    """
    imported = _lazy_import_plotting()
    if imported is None:
        return None
    np, plt = imported

    # (replicon, feature_type, mod_code) → summed raw_count
    agg: dict[tuple[str, str, str], int] = defaultdict(int)
    replicons: set[str] = set()
    feat_types: set[str] = set()
    for r in rows:
        agg[(r.feature.seqid, r.feature.feat_type, r.mod_code)] += r.raw_count
        replicons.add(r.feature.seqid)
        feat_types.add(r.feature.feat_type)

    if not agg:
        logger.warning("Per-feature plot: no methylated features — skipped.")
        return None

    replicons_s = sorted(replicons)
    feats_s     = sorted(feat_types)
    mods        = sorted(mod_codes)
    colours     = _mod_color_map(mod_codes, plt)

    nrows = len(replicons_s)
    fig, axes = plt.subplots(
        nrows, 1,
        figsize=(max(1.4 * len(feats_s) * max(len(mods), 1), 7), 3.0 * nrows),
        squeeze=False,
        sharex=True,
    )

    x        = np.arange(len(feats_s))
    n_mods   = max(len(mods), 1)
    width    = 0.8 / n_mods

    for i, rep in enumerate(replicons_s):
        ax = axes[i][0]
        for k, mod in enumerate(mods):
            heights = [agg.get((rep, ft, mod), 0) for ft in feats_s]
            offset  = (k - (n_mods - 1) / 2.0) * width
            label   = MOD_CODE_TO_LABEL.get(mod, mod)
            ax.bar(x + offset, heights, width=width,
                   color=colours[mod], label=f"{label} ({mod})")
        ax.set_ylabel("methylated\npositions", fontsize=8)
        ax.set_title(rep, fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels(feats_s, rotation=30, ha="right", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.margins(x=0.02)
        if i == 0:
            ax.legend(fontsize=8, title="modification", title_fontsize=8)

    fig.suptitle(
        "Methylation distribution per replicon and per feature, "
        "split by methylation type",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    fpath = output_dir / f"{stem}_per_replicon_feature.{fmt}"
    fig.savefig(fpath, dpi=dpi)
    plt.close(fig)
    logger.info("  [plot] per-replicon × feature distribution  →  %s", fpath.name)
    return fpath


def make_plots(
    features:          list[GffFeature],
    sites:             list[MethylSite],
    rows:              list[OutputRow],
    mod_codes:         set[str],
    declared_lengths:  dict[str, int],
    output_dir:        Path,
    output_name:       str,
    fmt:               str,
    bins:              int,
    dpi:               int,
) -> None:
    """
    Generate the figures:
        1. Always: methylation distribution along the replicon(s), split by
           methylation type.
        2. Only when more than one feature type carries methylation:
           the same distribution resolved per replicon and per feature type.
    """
    stem, _ext = _output_stem_and_ext(output_name)
    replicon_lengths = compute_replicon_lengths(features, sites, declared_lengths)

    logger.info("Generating figure(s):")
    plot_along_replicons(
        sites, replicon_lengths, mod_codes,
        output_dir, stem, fmt, bins, dpi,
    )

    distinct_feat_types = {r.feature.feat_type for r in rows}
    if len(distinct_feat_types) > 1:
        plot_per_replicon_per_feature(
            rows, mod_codes, output_dir, stem, fmt, dpi,
        )
    else:
        logger.info(
            "  [plot] per-feature plot skipped — only %d feature type with "
            "methylation.", len(distinct_feat_types),
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    gff_path   = Path(args.gff).resolve()
    output_dir = Path(args.output_dir).resolve()

    if not gff_path.is_file():
        logger.error("Annotation GFF file not found: %s", gff_path)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    output_name   = args.name or f"{gff_path.stem}_methylation_by_feature.tsv"
    feature_types = set(args.features) if args.features else None
    split         = not args.no_split

    # Determine methylation input source
    if args.bed:
        methyl_path   = Path(args.bed).resolve()
        methyl_source = f"BED  ({methyl_path})"
        loader        = load_bed
    else:
        methyl_path   = Path(args.methyl_gff).resolve()
        methyl_source = f"GFF  ({methyl_path})"
        loader        = load_methyl_gff

    if not methyl_path.is_file():
        logger.error("Methylation input file not found: %s", methyl_path)
        sys.exit(1)

    # Load explicit seqid map if provided
    explicit_map: dict[str, str] | None = None
    if args.seqid_map:
        explicit_map = {}
        with open(args.seqid_map, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    explicit_map[parts[0]] = parts[1]
        logger.info("Loaded %d explicit seqid mapping(s) from %s",
                    len(explicit_map), args.seqid_map)

    logger.info("Annotation GFF   : %s", gff_path)
    logger.info("Methylation input: %s", methyl_source)
    logger.info("Output directory : %s", output_dir)
    logger.info("Output base name : %s", output_name)
    logger.info("Split by mod type: %s", split)
    logger.info(
        "Feature types    : %s",
        ", ".join(sorted(feature_types)) if feature_types
        else f"all (excluding {', '.join(sorted(GFF_META_TYPES))})",
    )
    logger.info("Strand-aware     : %s", not args.ignore_strand)

    # Load data
    features = load_annotation_gff(gff_path, feature_types=feature_types)
    if not features:
        logger.error("No features loaded from annotation GFF. "
                     "Check --features or file content.")
        sys.exit(1)

    sites = loader(methyl_path)
    if not sites:
        logger.error("No methylation sites loaded. Check the input file.")
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Seqid reconciliation  ← KEY FIX
    # -----------------------------------------------------------------------
    if explicit_map is not None:
        seqid_map = explicit_map
    else:
        seqid_map = build_seqid_map(features, sites)

    sites = apply_seqid_map(sites, seqid_map)

    # Sanity check: warn if still no overlap after mapping
    gff_seqids = {f.seqid for f in features}
    bed_seqids = {s.chrom for s in sites}
    overlap    = gff_seqids & bed_seqids
    if not overlap:
        logger.error(
            "No chromosome names overlap between annotation (%s …) "
            "and methylation (%s …) after mapping. "
            "Supply a correct --seqid-map TSV (bed_name<TAB>gff_name).",
            ", ".join(sorted(gff_seqids)[:3]),
            ", ".join(sorted(bed_seqids)[:3]),
        )
        sys.exit(1)
    else:
        logger.info(
            "Seqid overlap after mapping: %d replicon(s) in common.", len(overlap)
        )

    # Index and cross-reference
    mod_codes    = {s.mod for s in sites}
    methyl_index = build_methyl_index(sites)

    logger.info(
        "Cross-referencing %d feature(s) × %d mod type(s) with %d site(s)…",
        len(features), len(mod_codes), len(sites),
    )
    rows = cross_reference(features, methyl_index, mod_codes, args.ignore_strand)

    # Write output
    logger.info("Writing output file(s):")
    write_outputs(rows, mod_codes, output_dir, output_name, split)

    # Figures
    if not args.no_plot:
        declared_lengths = load_sequence_regions(gff_path)
        make_plots(
            features, sites, rows, mod_codes, declared_lengths,
            output_dir, output_name, args.plot_format, args.bins, args.plot_dpi,
        )
    else:
        logger.info("Plotting disabled (--no-plot).")

    # Summary
    print_summary(features, rows, mod_codes, methyl_source)


if __name__ == "__main__":
    main()
