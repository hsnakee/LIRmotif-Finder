"""
Command-line interface for LIR Motif Finder.

Usage examples
--------------
# Basic run — writes results.csv + summary.txt to ./results/
lirmotiffinder input.fasta

# Custom thresholds and output directory
lirmotiffinder input.fasta --threshold 0.6 --min-disorder 10 --outdir my_results

# Only search for basic LIR motifs
lirmotiffinder input.fasta --motif-types basic

# Search basic and extended only
lirmotiffinder input.fasta --motif-types basic extended

# Output in TSV format
lirmotiffinder input.fasta --format tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
import textwrap
from pathlib import Path

from .core import run_analysis, ProteinResult, MotifHit, _MOTIF_PATTERNS


_VALID_MOTIF_TYPES = list(_MOTIF_PATTERNS.keys())


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
def write_csv(
    results: list[ProteinResult],
    out_path: Path,
    delimiter: str = ",",
) -> int:
    """Write all hits to a delimited file. Returns number of hits written."""
    fieldnames = [
        "protein_id",
        "protein_description",
        "protein_length",
        "motif_type",
        "motif_sequence",
        "motif_start_protein",
        "motif_end_protein",
        "disorder_start_protein",
        "disorder_end_protein",
        "disorder_length",
        "mean_disorder_score",
    ]
    hits_written = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for result in results:
            for hit in result.motif_hits:
                writer.writerow(hit.to_dict())
                hits_written += 1
    return hits_written


def write_summary(
    results: list[ProteinResult],
    out_path: Path,
    params: dict,
) -> None:
    """Write a plain-text summary."""
    total_proteins = len(results)
    proteins_with_hits = sum(1 for r in results if r.has_hits)
    total_hits = sum(len(r.motif_hits) for r in results)

    # Per-motif type counts
    type_counts: dict[str, int] = {}
    for result in results:
        for hit in result.motif_hits:
            type_counts[hit.motif_type] = type_counts.get(hit.motif_type, 0) + 1

    lines = [
        "=" * 60,
        "  LIR Motif Finder — Summary Report",
        "=" * 60,
        "",
        "Parameters",
        "----------",
    ]
    for k, v in params.items():
        lines.append(f"  {k}: {v}")
    lines += [
        "",
        "Results",
        "-------",
        f"  Total proteins analysed   : {total_proteins}",
        f"  Proteins with ≥1 LIR hit  : {proteins_with_hits}",
        f"  Total LIR hits            : {total_hits}",
        "",
        "Hits by motif type",
        "------------------",
    ]
    for mtype, count in sorted(type_counts.items()):
        lines.append(f"  {mtype:<20}: {count}")

    if proteins_with_hits:
        lines += [
            "",
            "Top 10 proteins by hit count",
            "----------------------------",
        ]
        ranked = sorted(results, key=lambda r: len(r.motif_hits), reverse=True)
        for r in ranked[:10]:
            if not r.has_hits:
                break
            lines.append(f"  {r.protein_id:<20} {len(r.motif_hits)} hit(s)")

    lines += ["", "=" * 60]

    out_path.write_text("\n".join(lines))


def print_hits_to_stdout(results: list[ProteinResult]) -> None:
    """Print a brief human-readable table to stdout."""
    header = (
        f"{'Protein ID':<20} {'Motif type':<18} {'Motif':<8} "
        f"{'Position (protein)':<22} {'Disorder region':<20} {'Score':>6}"
    )
    print(header)
    print("-" * len(header))
    for result in results:
        for hit in result.motif_hits:
            pos = f"{hit.motif_start_protein}-{hit.motif_end_protein}"
            dis = f"{hit.disorder_start_protein}-{hit.disorder_end_protein}"
            print(
                f"{hit.protein_id:<20} {hit.motif_type:<18} {hit.motif_sequence:<8} "
                f"{pos:<22} {dis:<20} {hit.mean_disorder_score:>6.3f}"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lirmotiffinder",
        description=textwrap.dedent(
            """\
            LIR Motif Finder
            ================
            Identifies LIR (LC3-Interacting Region) motifs within
            intrinsically disordered regions predicted by metapredict.

            Motif patterns searched:
              basic           [WFY]xx[LIV]
              extended        [DE][WFY]xx[LIV]
              acidic_extended [DE]{2}[WFY]xx[LIV]
            """
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "fasta",
        metavar="FASTA",
        help="Input FASTA file (protein sequences).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        metavar="FLOAT",
        help="Metapredict disorder score threshold (default: 0.5).",
    )
    parser.add_argument(
        "--min-disorder",
        type=int,
        default=5,
        metavar="INT",
        help="Minimum disordered region length in residues (default: 5).",
    )
    parser.add_argument(
        "--motif-types",
        nargs="+",
        choices=_VALID_MOTIF_TYPES,
        default=None,
        metavar="TYPE",
        help=(
            f"Motif types to search for. Choose from: "
            f"{', '.join(_VALID_MOTIF_TYPES)}. Default: all."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results"),
        metavar="DIR",
        help="Output directory (default: ./results/).",
    )
    parser.add_argument(
        "--format",
        choices=["csv", "tsv"],
        default="csv",
        help="Output format for the hits table (default: csv).",
    )
    parser.add_argument(
        "--no-summary",
        action="store_true",
        help="Skip writing the summary text file.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-hit stdout output.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        parser.error(f"FASTA file not found: {fasta_path}")

    # Set up output directory
    args.outdir.mkdir(parents=True, exist_ok=True)
    delimiter = "," if args.format == "csv" else "\t"
    hits_filename = f"lir_hits.{args.format}"

    params = {
        "fasta": str(fasta_path),
        "disorder_threshold": args.threshold,
        "min_disorder_length": args.min_disorder,
        "motif_types": args.motif_types or "all",
        "output_directory": str(args.outdir),
    }

    print(f"[lirmotiffinder] Analysing: {fasta_path}", file=sys.stderr)
    print(f"[lirmotiffinder] Disorder threshold: {args.threshold}", file=sys.stderr)
    print(
        f"[lirmotiffinder] Min disorder region length: {args.min_disorder} residues",
        file=sys.stderr,
    )

    results = run_analysis(
        fasta_path=str(fasta_path),
        threshold=args.threshold,
        min_disorder_length=args.min_disorder,
        motif_types=args.motif_types,
    )

    total_hits = sum(len(r.motif_hits) for r in results)
    print(
        f"[lirmotiffinder] Done. {len(results)} proteins, {total_hits} LIR hit(s).",
        file=sys.stderr,
    )

    # Write hits table
    hits_path = args.outdir / hits_filename
    n = write_csv(results, hits_path, delimiter=delimiter)
    print(f"[lirmotiffinder] Hits written to: {hits_path} ({n} rows)", file=sys.stderr)

    # Write summary
    if not args.no_summary:
        summary_path = args.outdir / "summary.txt"
        write_summary(results, summary_path, params)
        print(f"[lirmotiffinder] Summary written to: {summary_path}", file=sys.stderr)

    # Print to stdout
    if not args.quiet and total_hits > 0:
        print()
        print_hits_to_stdout(results)


if __name__ == "__main__":
    main()
