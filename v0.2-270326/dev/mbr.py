#!/usr/bin/env python3

"""
mbr.py
-------------------
Filters a bedMethyl file produced by modkit (ONT)
and saves a BED file containing only the first 6 columns:
    1. chrom        - chromosome / contig name
    2. start        - 0-based start position
    3. end          - 0-based end position (start + 1 for CpG)
    4. modified_base_code_and_motif - e.g. "m" for 5mC
    5. score        - coverage score (0-1000)
    6. strand       - "+" or "-"

Features
--------
  - Custom output file name  (-n / --name)
  - Optional split by replicon (--split-by-replicon)
  - Optional quality / score filtering (--min-score)
  - Summary printed to stdout:
      replicons found, methylation types and their occurrence counts

Usage
-----
    # Basic
    python3 mbr.py -i input.bed -o results/

    # Custom output name
    python3 mbr.py -i input.bed -o results/ -n my_output.bed

    # Split into one file per replicon
    python3 mbr.py -i input.bed -o results/ --split-by-replicon

    # Keep only high-confidence sites (score >= 500)
    python3 mbr.py -i input.bed -o results/ --min-score 500

    # Combine all options
    python3 mbr.py -i input.bed -o results/ --split-by-replicon --min-score 300

"""

import argparse
import sys
import logging
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Column indices (0-based) in a modkit bedMethyl file
# ---------------------------------------------------------------------------
CHROM_COL       = 0    # chromosome / replicon
START_COL       = 1    # start position
END_COL         = 2    # end position
MOD_BASE_COL    = 3    # modified base code & motif  (e.g. "m", "h", "a")
SCORE_COL       = 4    # score (0–1000)
STRAND_COL      = 5    # strand (+ / -)

N_REQUIRED_COLS = 6    # minimum columns expected per data line

# Human-readable labels for common modkit modification codes
MOD_BASE_LABELS = {
    "m":  "5-methylcytosine (5mC)",
    "h":  "5-hydroxymethylcytosine (5hmC)",
    "a":  "6-methyladenine (6mA)",
    "f":  "5-formylcytosine (5fC)",
    "c":  "5-carboxylcytosine (5caC)",
    "g":  "N7-methylguanine (7mG)",
    "e":  "N3-methylcytosine (3mC)",
    "b":  "5-hydroxymethyluracil (5hmU)",
}


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------
def print_summary(
    replicon_counts: dict,
    mod_counts: dict,
    records_written: int,
    records_skipped: int,
    records_filtered: int,
) -> None:
    """Print a formatted summary table to stdout."""
    sep  = "=" * 62
    sep2 = "-" * 62

    print(f"\n{sep}")
    print(" Summary")
    print(sep)

    # --- Replicons ---
    print(f"\n  Replicons found: {len(replicon_counts)}")
    print(f"  {sep2}")
    print(f"  {'Replicon':<35} {'Records':>10}")
    print(f"  {sep2}")
    for replicon, count in sorted(replicon_counts.items()):
        print(f"  {replicon:<35} {count:>10,}")

    # --- Methylation types ---
    print(f"\n  Methylation types found: {len(mod_counts)}")
    print(f"  {sep2}")
    print(f"  {'Type':<6}  {'Description':<38} {'Occurrences':>10}")
    print(f"  {sep2}")
    for mod, count in sorted(mod_counts.items(), key=lambda x: -x[1]):
        label = MOD_BASE_LABELS.get(mod, "unknown modification")
        print(f"  {mod:<6}  {label:<38} {count:>10,}")

    # --- Record counts ---
    print(f"\n  Record statistics:")
    print(f"  {sep2}")
    print(f"  {'Written':<30} {records_written:>10,}")
    print(f"  {'Filtered out (score threshold)':<30} {records_filtered:>10,}")
    print(f"  {'Skipped (malformed lines)':<30} {records_skipped:>10,}")
    print(f"\n{sep}\n")


# ---------------------------------------------------------------------------
# Core filtering function
# ---------------------------------------------------------------------------
def filter_bedmethyl(
    input_path: Path,
    output_dir: Path,
    output_name: str | None = None,
    split_by_replicon: bool = False,
    min_score: int | None = None,
    keep_headers: bool = True,
    comment_char: str = "#",
) -> None:
    """
    Parameters
    ----------
    input_path : Path
        Input bedMethyl file.
    output_dir : Path
        Destination directory (must already exist).
    output_name : str | None
        Base filename for single-file output. Defaults to
        '<input_stem>_filtered.bed'. Ignored when split_by_replicon=True.
    split_by_replicon : bool
        Write one BED file per chromosome/contig.
    min_score : int | None
        Minimum score (0–1000). Records below this threshold are discarded.
    keep_headers : bool
        Copy comment/header lines to all output files.
    comment_char : str
        Character that marks comment/header lines.
    """

    stem = input_path.stem

    header_lines: list[str] = []

    # Per-replicon output lines and statistics
    replicon_data: dict[str, list[str]] = defaultdict(list)
    replicon_counts: dict[str, int]     = defaultdict(int)
    mod_counts: dict[str, int]          = defaultdict(int)

    records_written  = 0
    records_skipped  = 0
    records_filtered = 0

    # -----------------------------------------------------------------------
    # Read and filter
    # -----------------------------------------------------------------------
    with open(input_path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")

            if not raw:
                continue

            # Header / comment lines
            if raw.startswith(comment_char):
                header_lines.append(raw)
                continue

            cols = raw.split("\t")

            # Column count validation
            if len(cols) < N_REQUIRED_COLS:
                logger.warning(
                    "Line %d: only %d column(s) found (expected >= %d) — skipped.",
                    line_no, len(cols), N_REQUIRED_COLS,
                )
                records_skipped += 1
                continue

            chrom    = cols[CHROM_COL]
            mod_base = cols[MOD_BASE_COL]

            # Parse and validate score
            try:
                score = int(cols[SCORE_COL])
            except ValueError:
                logger.warning(
                    "Line %d: cannot parse score '%s' as integer — skipped.",
                    line_no, cols[SCORE_COL],
                )
                records_skipped += 1
                continue

            # Score threshold filter
            if min_score is not None and score < min_score:
                records_filtered += 1
                continue

            # Build 6-column output line
            out_line = "\t".join([
                chrom,
                cols[START_COL],
                cols[END_COL],
                mod_base,
                str(score),
                cols[STRAND_COL],
            ])

            replicon_data[chrom].append(out_line)
            replicon_counts[chrom] += 1
            mod_counts[mod_base]   += 1
            records_written        += 1

    # -----------------------------------------------------------------------
    # Write output
    # -----------------------------------------------------------------------
    if split_by_replicon:
        logger.info(
            "Split-by-replicon mode: writing %d file(s).", len(replicon_data)
        )
        for replicon, lines in sorted(replicon_data.items()):
            # Sanitise replicon name for safe use in filenames
            safe_name = (
                replicon.replace("/", "_")
                        .replace("\\", "_")
                        .replace(" ", "_")
            )
            fname = f"{stem}_{safe_name}_filtered.bed"
            fpath = output_dir / fname

            with open(fpath, "w", encoding="utf-8") as fh:
                if keep_headers:
                    for h in header_lines:
                        fh.write(h + "\n")
                for record in lines:
                    fh.write(record + "\n")

            logger.info("  %6d record(s)  →  %s", len(lines), fpath.name)

    else:
        # Single output file
        if output_name is None:
            output_name = f"{stem}_filtered.bed"

        fpath = output_dir / output_name

        with open(fpath, "w", encoding="utf-8") as fh:
            if keep_headers:
                for h in header_lines:
                    fh.write(h + "\n")
            for lines in replicon_data.values():
                for record in lines:
                    fh.write(record + "\n")

        logger.info("Filtered BED file saved to: %s", fpath)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print_summary(
        replicon_counts=replicon_counts,
        mod_counts=mod_counts,
        records_written=records_written,
        records_skipped=records_skipped,
        records_filtered=records_filtered,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a modkit bedMethyl file to retain only the first 6 BED columns\n"
            "(chrom, start, end, modified_base, score, strand).\n"
            "Optionally split output by replicon and/or apply a score threshold."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  Basic:\n"
            "    python3 mbr.py -i input.bed -o results/\n\n"
            "  Custom output name:\n"
            "    python3 mbr.py -i input.bed -o results/ -n my_output.bed\n\n"
            "  Split into one file per replicon:\n"
            "    python3 mbr.py -i input.bed -o results/ --split-by-replicon\n\n"
            "  Strict (score >= 500):\n"
            "    python3 mbr.py -i input.bed -o results/ --min-score 500\n\n"
            "  Wizard mode:\n"
            "    python3 mbr.py -i input.bed -o results/ \\\n"
            "        --split-by-replicon --min-score 300\n"
        ),
    )

    # Required arguments
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-i", "--input",
        required=True,
        metavar="FILE",
        help="Path to the input modkit bedMethyl file.",
    )
    req.add_argument(
        "-o", "--output-dir",
        required=True,
        metavar="DIR",
        help="Directory where output BED file(s) will be saved (created if absent).",
    )

    # Output options
    out = parser.add_argument_group("output options")
    out.add_argument(
        "-n", "--name",
        default=None,
        metavar="FILENAME",
        help=(
            "Name for the output BED file (e.g. 'results.bed'). "
            "Defaults to '<input_stem>_filtered.bed'. "
            "Ignored when --split-by-replicon is active."
        ),
    )
    out.add_argument(
        "--split-by-replicon",
        action="store_true",
        help=(
            "Write one BED file per replicon (chromosome/contig). "
            "Output files are named '<input_stem>_<replicon>_filtered.bed'."
        ),
    )
    out.add_argument(
        "--no-header-passthrough",
        action="store_true",
        help="Do NOT copy comment/header lines (starting with '#') to output files.",
    )

    # Filtering options
    flt = parser.add_argument_group("filtering options")
    flt.add_argument(
        "--min-score",
        type=int,
        default=None,
        metavar="INT",
        help=(
            "Minimum score threshold (0–1000). Records whose score is strictly "
            "below this value are excluded from all outputs. "
            "The score column in modkit bedMethyl represents coverage depth "
            "rescaled to the 0–1000 range."
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    input_path = Path(args.input).resolve()
    output_dir = Path(args.output_dir).resolve()

    # Validate input file
    if not input_path.is_file():
        logger.error("Input file not found: %s", input_path)
        sys.exit(1)

    # Validate score threshold
    if args.min_score is not None and not (0 <= args.min_score <= 1000):
        logger.error(
            "--min-score must be an integer between 0 and 1000 (got %d).",
            args.min_score,
        )
        sys.exit(1)

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Log parameters
    logger.info("Input            : %s", input_path)
    logger.info("Output directory : %s", output_dir)
    if args.name and not args.split_by_replicon:
        logger.info("Output name      : %s", args.name)
    if args.split_by_replicon:
        logger.info("Split by replicon: enabled")
    if args.min_score is not None:
        logger.info("Min score filter : >= %d", args.min_score)

    filter_bedmethyl(
        input_path=input_path,
        output_dir=output_dir,
        output_name=args.name,
        split_by_replicon=args.split_by_replicon,
        min_score=args.min_score,
        keep_headers=not args.no_header_passthrough,
    )


if __name__ == "__main__":
    main()