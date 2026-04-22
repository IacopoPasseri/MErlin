#!/usr/bin/env python3

"""
bin.py
----------------------
Counts methylation events falling inside fixed-size genomic bins (or a sliding
window with a user-defined step) across every replicon found in the input file.

Methylation input — choose ONE of:
    -b / --bed          Filtered 6-column BED file from xam.py
                        (ONT/modkit output)
    -m / --methyl-gff   GFF3 file produced from a PacBio sequencing basecalling,
                        where each feature record represents one methylation call.


Genome sizes — choose ONE of:
    -g / --genome-size  Two-column TSV  <replicon>  <length_bp>  (no header)
                        When provided, bins extend to the true chromosome end
                        rather than the last observed methylation site.

Windowing modes
---------------
    Fixed bins (default)
        The genome is partitioned into non-overlapping bins of --bin-width bp.
        Each bin [start, start + bin_width) is counted independently.

    Sliding window  (--step < --bin-width)
        Windows of --bin-width bp are placed every --step bp.
        Windows overlap; a methylation site may be counted in multiple windows.

Output
------
A tab-separated BED-like file with one row per (replicon, bin, strand, mod_type)
that contains at least one methylation hit:

    replicon  bin_start  bin_end  strand  mod_code  mod_label  raw_count  norm_count

By default one output file is written per modification type, named
    <stem>_<MOD_LABEL>.bed
Use --no-split to write a single combined file.


Strand handling
---------------
By default counts are reported separately for each strand (+, -) within every
bin.  Use --ignore-strand to merge both strands into a single row per bin.

Usage
-----
    # ONT BED input, 5 000 bp bins
    python3 bin.py -b filtered.bed -o results/ --bin-width 5000

    # PacBio GFF input, 10 000 bp bins, true genome sizes
    python3 bin.py -m basecalling.gff -o results/ \
        --bin-width 10000 --genome-size sizes.tsv

    # Sliding window: 5 000 bp window, 1 000 bp step
    python3 bin.py -b filtered.bed -o results/ \
        --bin-width 5000 --step 1000

    # Single output file, ignore strand
    python3 bin.py -b filtered.bed -o results/ \
        --bin-width 5000 --no-split --ignore-strand

    # Custom output name
    python3 bin.py -b filtered.bed -o results/ \
        --bin-width 5000 -n my_windows.bed

"""

import argparse
import bisect
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
# Constants  (identical to annotate_methylation_gff.py for consistency)
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

MOD_LABEL_TO_CODE: dict[str, str] = {v: k for k, v in MOD_CODE_TO_LABEL.items()}

METHYL_GFF_TYPE_SYNONYMS: dict[str, str] = {
    "cpg":           "m",
    "5mc":           "m",
    "5-mc":          "m",
    "5methylc":      "m",
    "methylation":   "m",
    "modified_base": "m",
    "6ma":           "a",
    "6-ma":          "a",
    "6methyla":      "a",
    "dam":           "a",
    "gatc":          "a",
    "dcm":           "m",
    "ccwgg":         "m",
    "5hmc":          "h",
    "5-hmc":         "h",
}

GFF_META_TYPES: set[str] = {"region", "sequence_feature", "sequence_alteration"}

# GFF3 column indices (0-based)
GFF_SEQID  = 0
GFF_TYPE   = 2
GFF_START  = 3   # 1-based inclusive
GFF_STRAND = 6
GFF_ATTRS  = 8

# BED column indices (0-based)
BED_CHROM  = 0
BED_START  = 1   # 0-based
BED_MOD    = 3
BED_STRAND = 5

OUTPUT_HEADER = "\t".join([
    "replicon",
    "bin_start",    # 0-based
    "bin_end",      # 0-based exclusive
    "strand",
    "mod_code",
    "mod_label",
    "raw_count",
    "norm_count",   # raw_count / bin_size
])


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MethylSite:
    """One methylation call from BED or basecalling GFF."""
    chrom:  str
    start:  int    # 0-based
    mod:    str    # normalised mod code
    strand: str


# ---------------------------------------------------------------------------
# Shared GFF3 helpers
# ---------------------------------------------------------------------------

## Strip the attribute column to get the key elements
def _parse_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            key, _, val = part.partition("=")
            attrs[key.strip()] = val.strip()
    return attrs

## Tab separation
def _split_gff_line(raw: str) -> list[str] | None:
    cols = raw.split("\t")
    if len(cols) < 9:
        cols = raw.split()
    return cols if len(cols) >= 9 else None


# ---------------------------------------------------------------------------
# Mod-code resolver  (for GFF input)
# ---------------------------------------------------------------------------

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
        if val.lower() in METHYL_GFF_TYPE_SYNONYMS:
            return METHYL_GFF_TYPE_SYNONYMS[val.lower()]
    return feat_type


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

## Load methylation sites from a 6-column filtered BED file (ONT/modkit)
def load_bed(bed_path: Path) -> list[MethylSite]:
    """
    Columns used (0-based index):
        0  chrom
        1  start  (0-based)
        3  mod code
        5  strand
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
        "BED: %d site(s) loaded  |  %d skipped (malformed).",
        len(sites), skipped,
    )
    return sites

## Load methylation sites from a PacBio basecalling GFF3 file
def load_methyl_gff(gff_path: Path) -> list[MethylSite]:
    """
    Each data line represents one methylation call:
        col 0  seqid   — replicon
        col 2  type    — mod type (resolved via _resolve_mod_code)
        col 3  start   — 1-based inclusive -> converted to 0-based
        col 6  strand
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
            if feat_type in GFF_META_TYPES:
                skipped_meta += 1
                continue

            attrs    = _parse_attributes(cols[GFF_ATTRS])
            mod_code = _resolve_mod_code(feat_type, attrs)

            if mod_code not in MOD_CODE_TO_LABEL and mod_code not in unknown_types:
                unknown_types.add(mod_code)
                logger.warning(
                    "Unrecognised modification type '%s' — kept as-is.", mod_code
                )

            try:
                start_0 = int(cols[GFF_START]) - 1   # GFF 1-based -> 0-based
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
        "GFF: %d site(s) loaded  |  %d skipped (malformed)  |  %d skipped (meta).",
        len(sites), skipped_bad, skipped_meta,
    )
    return sites

## Load replicon sizes from a two-column TSV (no header)
def load_genome_sizes(path: Path) -> dict[str, int]:
    """
        <replicon_name>  <length_bp>
    """
    sizes: dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split("\t")
            if len(parts) < 2:
                parts = raw.split()
            if len(parts) < 2:
                logger.warning(
                    "Genome-size line %d: expected 2 columns — skipped.", line_no
                )
                continue
            try:
                sizes[parts[0]] = int(parts[1])
            except ValueError:
                logger.warning(
                    "Genome-size line %d: cannot parse length '%s' — skipped.",
                    line_no, parts[1],
                )
    logger.info("Genome sizes loaded for %d replicon(s).", len(sizes))
    return sizes


# ---------------------------------------------------------------------------
# Position index for O(log n) range queries
# ---------------------------------------------------------------------------

## chrom, strand, mod_code -> sorted list of 0-based positions
def build_index(
    sites: list[MethylSite],
) -> dict[tuple[str, str, str], list[int]]:
    
    index: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    for site in sites:
        index[(site.chrom, site.strand, site.mod)].append(site.start)
    for key in index:
        index[key].sort()
    return index

## Count sorted positions in [start, end) using binary search
def count_in_range(positions: list[int], start: int, end: int) -> int:
    lo = bisect.bisect_left(positions, start)
    hi = bisect.bisect_left(positions, end, lo)
    return hi - lo


# ---------------------------------------------------------------------------
# Output row
# ---------------------------------------------------------------------------

@dataclass
class BinRow:
    replicon:  str
    bin_start: int    # 0-based
    bin_end:   int    # 0-based exclusive
    strand:    str
    mod_code:  str
    raw_count: int

    @property
    def norm_count(self) -> float:
        size = self.bin_end - self.bin_start
        return self.raw_count / size if size > 0 else 0.0

    def to_bed(self) -> str:
        label = MOD_CODE_TO_LABEL.get(self.mod_code, self.mod_code)
        return "\t".join([
            self.replicon,
            str(self.bin_start),
            str(self.bin_end),
            self.strand,
            self.mod_code,
            label,
            str(self.raw_count),
            f"{self.norm_count:.8f}",
        ])


# ---------------------------------------------------------------------------
# Core binning
# ---------------------------------------------------------------------------

def bin_methylation(
    sites:         list[MethylSite],
    genome_sizes:  dict[str, int] | None,
    bin_width:     int,
    step:          int,
    ignore_strand: bool,
) -> list[BinRow]:
    """
    Partition every replicon into windows and count methylation per bin.

    Parameters
    ----------
    sites         : loaded methylation sites.
    genome_sizes  : optional dict replicon -> length_bp.
                    If None, lengths are inferred from max(site.start) + 1.
    bin_width     : window size in bp.
    step          : window step in bp (== bin_width -> non-overlapping bins).
    ignore_strand : if True, '+' and '-' sites are merged into a single '.' row.

    Returns
    -------
    list[BinRow]  — only bins with raw_count > 0 are included.
    """
    index     = build_index(sites)
    mod_codes = sorted({s.mod for s in sites})

    # Build replicon length map ---------------------------------------------------
    # Start from provided genome sizes (if any), then extend/correct from data.
    replicons: dict[str, int] = {}
    if genome_sizes:
        replicons.update(genome_sizes)

    # Ensure every replicon that has sites is represented and length is >= data max
    for site in sites:
        chrom    = site.chrom
        observed = site.start + 1
        if chrom not in replicons:
            replicons[chrom] = observed
        elif replicons[chrom] < observed:
            # Data extends beyond the provided genome size — warn and extend
            logger.warning(
                "Replicon '%s': site at position %d exceeds provided genome "
                "size %d bp — extending to cover the data.",
                chrom, site.start, replicons[chrom],
            )
            replicons[chrom] = observed

    if not replicons:
        logger.error("No replicons found in input data.")
        return []

    # Strands to query individually (or merged)
    query_strands = ["+", "-", "."]   # always query all three index keys

    rows: list[BinRow] = []
    total_bins = 0

    for chrom in sorted(replicons):
        chrom_len = replicons[chrom]

        if chrom_len == 0:
            continue

        # Generate bin start positions
        bin_starts = range(0, chrom_len, step)

        for bin_start in bin_starts:
            bin_end = min(bin_start + bin_width, chrom_len)
            total_bins += 1

            for mod in mod_codes:
                if ignore_strand:
                    # Aggregate all strands into a single '.' row
                    total = 0
                    for strand in query_strands:
                        positions = index.get((chrom, strand, mod))
                        if positions:
                            total += count_in_range(positions, bin_start, bin_end)
                    if total > 0:
                        rows.append(BinRow(
                            replicon  = chrom,
                            bin_start = bin_start,
                            bin_end   = bin_end,
                            strand    = ".",
                            mod_code  = mod,
                            raw_count = total,
                        ))
                else:
                    # Separate rows for each strand that has hits
                    for strand in query_strands:
                        positions = index.get((chrom, strand, mod))
                        if not positions:
                            continue
                        count = count_in_range(positions, bin_start, bin_end)
                        if count > 0:
                            rows.append(BinRow(
                                replicon  = chrom,
                                bin_start = bin_start,
                                bin_end   = bin_end,
                                strand    = strand,
                                mod_code  = mod,
                                raw_count = count,
                            ))

    logger.info(
        "Binning: %d replicon(s)  |  %d bin(s) evaluated  |  "
        "%d (bin x strand x mod) rows with hits",
        len(replicons), total_bins, len(rows),
    )
    return rows


# ---------------------------------------------------------------------------
# Output writer
# ---------------------------------------------------------------------------

def _stem_and_ext(output_name: str) -> tuple[str, str]:
    p   = Path(output_name)
    ext = p.suffix if p.suffix else ".bed"
    return p.stem, ext


def write_outputs(
    rows:        list[BinRow],
    mod_codes:   set[str],
    output_dir:  Path,
    output_name: str,
    split:       bool,
) -> None:
    """
    split=True  (default): one file per mod type -> <stem>_<MOD_LABEL><ext>
    split=False           : single combined file -> <output_name>
    """
    stem, ext = _stem_and_ext(output_name)

    if split:
        for mod in sorted(mod_codes):
            label    = MOD_CODE_TO_LABEL.get(mod, mod)
            fpath    = output_dir / f"{stem}_{label}{ext}"
            mod_rows = [r for r in rows if r.mod_code == mod]
            with open(fpath, "w", encoding="utf-8") as fh:
                fh.write("# " + OUTPUT_HEADER + "\n")
                for row in mod_rows:
                    fh.write(row.to_bed() + "\n")
            logger.info(
                "  [%s]  %d row(s)  ->  %s", label, len(mod_rows), fpath.name
            )
    else:
        fpath = output_dir / output_name
        with open(fpath, "w", encoding="utf-8") as fh:
            fh.write("# " + OUTPUT_HEADER + "\n")
            for row in rows:
                fh.write(row.to_bed() + "\n")
        logger.info("Combined output: %s  (%d rows)", fpath, len(rows))


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------

def print_summary(
    sites:     list[MethylSite],
    rows:      list[BinRow],
    bin_width: int,
    step:      int,
) -> None:
    sep  = "=" * 64
    sep2 = "-" * 64

    replicon_sites: dict[str, int] = defaultdict(int)
    replicon_hits:  dict[str, int] = defaultdict(int)
    mod_sites:      dict[str, int] = defaultdict(int)
    mod_hits:       dict[str, int] = defaultdict(int)

    for s in sites:
        replicon_sites[s.chrom] += 1
        mod_sites[s.mod]        += 1

    for r in rows:
        replicon_hits[r.replicon] += r.raw_count
        mod_hits[r.mod_code]      += r.raw_count

    mode_str = (
        "fixed non-overlapping bins"
        if step == bin_width
        else f"sliding window  (step = {step:,} bp)"
    )

    print(f"\n{sep}")
    print("  Methylation Windowing — Summary")
    print(sep)
    print(f"\n  Bin width    : {bin_width:,} bp")
    print(f"  Mode         : {mode_str}")
    print(f"  Output rows  : {len(rows):,}  (bins with >= 1 methylation hit)")

    print(f"\n  Replicons:")
    print(f"  {sep2}")
    print(f"  {'Replicon':<36} {'Total sites':>12} {'Sum hits in bins':>17}")
    print(f"  {sep2}")
    for rep in sorted(replicon_sites):
        print(
            f"  {rep:<36} {replicon_sites[rep]:>12,} "
            f"{replicon_hits.get(rep, 0):>17,}"
        )

    print(f"\n  Modification types:")
    print(f"  {sep2}")
    print(f"  {'Code':<6}  {'Label':<8}  {'Total sites':>12} {'Sum hits':>10}")
    print(f"  {sep2}")
    for mod in sorted(mod_sites, key=lambda m: -mod_sites[m]):
        label = MOD_CODE_TO_LABEL.get(mod, mod)
        print(
            f"  {mod:<6}  {label:<8}  {mod_sites[mod]:>12,} "
            f"{mod_hits.get(mod, 0):>10,}"
        )

    print(f"\n{sep}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count methylation events in fixed bins or sliding windows.\n"
            "Accepts ONT BED or PacBio GFF input.\n\n"
            "Provide exactly one of -b/--bed or -m/--methyl-gff."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  ONT BED, 5 000 bp fixed bins:\n"
            "    python bin.py \\\n"
            "        -b filtered.bed -o results/ --bin-width 5000\n\n"
            "  PacBio GFF, 10 000 bp bins, true genome sizes:\n"
            "    python bin.py \\\n"
            "        -m basecalling.gff -o results/ \\\n"
            "        --bin-width 10000 --genome-size sizes.tsv\n\n"
            "  Sliding window (5 000 bp window, 1 000 bp step):\n"
            "    python bin.py \\\n"
            "        -b filtered.bed -o results/ \\\n"
            "        --bin-width 5000 --step 1000\n\n"
            "  Single combined output, both strands merged:\n"
            "    python bin.py \\\n"
            "        -b filtered.bed -o results/ \\\n"
            "        --bin-width 5000 --no-split --ignore-strand\n"
        ),
    )

    # Methylation input (mutually exclusive, one required)
    meth    = parser.add_argument_group("methylation input (provide exactly one)")
    meth_ex = meth.add_mutually_exclusive_group(required=True)
    meth_ex.add_argument(
        "-b", "--bed", metavar="FILE",
        help="Filtered 6-column ONT BED file from filter_bedmethyl.py.",
    )
    meth_ex.add_argument(
        "-m", "--methyl-gff", metavar="FILE",
        help="PacBio basecalling GFF3 file (one methylation call per line).",
    )

    # Required output
    req = parser.add_argument_group("required output")
    req.add_argument(
        "-o", "--output-dir", required=True, metavar="DIR",
        help="Directory where output BED file(s) will be saved.",
    )

    # Windowing
    win = parser.add_argument_group("windowing options")
    win.add_argument(
        "--bin-width", type=int, required=True, metavar="INT",
        help="Width of each genomic window / bin in base pairs.",
    )
    win.add_argument(
        "--step", type=int, default=None, metavar="INT",
        help=(
            "Step between consecutive window starts (bp). "
            "Defaults to --bin-width (non-overlapping fixed bins). "
            "Use a value smaller than --bin-width for a sliding window."
        ),
    )

    # Genome sizes
    parser.add_argument(
        "-g", "--genome-size", default=None, metavar="FILE",
        help=(
            "Two-column TSV: <replicon>  <length_bp>  (no header). "
            "If omitted, lengths are inferred from the last observed "
            "methylation position + 1."
        ),
    )

    # Output options
    out = parser.add_argument_group("output options")
    out.add_argument(
        "-n", "--name", default=None, metavar="FILENAME",
        help=(
            "Base name for the output file(s). "
            "Defaults to '<input_stem>_w<bin_width>.bed'. "
            "When splitting (default), the mod label is appended before the extension."
        ),
    )
    out.add_argument(
        "--no-split", action="store_true",
        help="Write a single combined BED instead of one file per mod type.",
    )

    # Strand
    flt = parser.add_argument_group("strand options")
    flt.add_argument(
        "--ignore-strand", action="store_true",
        help=(
            "Merge + and - strand sites into a single count per bin (strand = '.'). "
            "Default: separate rows for each strand."
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir).resolve()

    # Resolve methylation input
    if args.bed:
        methyl_path  = Path(args.bed).resolve()
        loader       = load_bed
        source_label = f"BED  ({methyl_path.name})"
    else:
        methyl_path  = Path(args.methyl_gff).resolve()
        loader       = load_methyl_gff
        source_label = f"GFF  ({methyl_path.name})"

    if not methyl_path.is_file():
        logger.error("Methylation input file not found: %s", methyl_path)
        sys.exit(1)

    # Validate windowing parameters
    bin_width = args.bin_width
    step      = args.step if args.step is not None else bin_width

    if bin_width <= 0:
        logger.error("--bin-width must be a positive integer.")
        sys.exit(1)
    if step <= 0:
        logger.error("--step must be a positive integer.")
        sys.exit(1)
    if step > bin_width:
        logger.error(
            "--step (%d) cannot exceed --bin-width (%d).", step, bin_width
        )
        sys.exit(1)

    # Default output name
    output_name = args.name or f"{methyl_path.stem}_w{bin_width}.bed"
    split       = not args.no_split

    output_dir.mkdir(parents=True, exist_ok=True)

    # Genome sizes (optional)
    genome_sizes: dict[str, int] | None = None
    if args.genome_size:
        gs_path = Path(args.genome_size).resolve()
        if not gs_path.is_file():
            logger.error("Genome-size file not found: %s", gs_path)
            sys.exit(1)
        genome_sizes = load_genome_sizes(gs_path)

    # Log run parameters
    logger.info("Methylation input : %s", source_label)
    logger.info("Output directory  : %s", output_dir)
    logger.info("Output base name  : %s", output_name)
    logger.info("Bin width         : %d bp", bin_width)
    if step == bin_width:
        logger.info("Mode              : fixed non-overlapping bins")
    else:
        logger.info("Mode              : sliding window  (step = %d bp)", step)
    logger.info("Split by mod type : %s", split)
    logger.info("Strand-aware      : %s", not args.ignore_strand)
    if genome_sizes:
        logger.info("Genome sizes      : %d replicon(s) from file", len(genome_sizes))
    else:
        logger.info("Genome sizes      : inferred from data  (max position + 1)")

    # Load methylation sites
    sites = loader(methyl_path)
    if not sites:
        logger.error("No methylation sites loaded. Check input file.")
        sys.exit(1)

    # Bin
    n_replicons = len({s.chrom for s in sites})
    logger.info(
        "Binning %d site(s) across %d replicon(s)...",
        len(sites), n_replicons,
    )
    rows = bin_methylation(
        sites         = sites,
        genome_sizes  = genome_sizes,
        bin_width     = bin_width,
        step          = step,
        ignore_strand = args.ignore_strand,
    )

    if not rows:
        logger.warning(
            "No bins contained any methylation sites. "
            "Check --bin-width relative to data density."
        )
        sys.exit(0)

    # Write
    mod_codes = {r.mod_code for r in rows}
    logger.info("Writing output file(s):")
    write_outputs(rows, mod_codes, output_dir, output_name, split)

    # Summary
    print_summary(sites, rows, bin_width, step)


if __name__ == "__main__":
    main()