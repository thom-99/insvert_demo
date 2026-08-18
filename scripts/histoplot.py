#!/usr/bin/env python3
"""Plot supported structural-variant record counts per chromosome."""

import argparse
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam


SV_TYPES = ("INS", "DEL", "INV", "DUP", "BND")
SV_COLORS = {
    "INS": "#4C78A8",
    "DEL": "#F58518",
    "INV": "#B279A2",
    "DUP": "#54A24B",
    "BND": "#E45756",
}


def natural_chromosome_key(chromosome):
    """Return a key that orders chr2 before chr10 and numeric contigs before text."""
    label = chromosome.casefold()
    if label.startswith("chr"):
        label = label[3:]

    return tuple(
        (0, int(part)) if part.isdigit() else (1, part)
        for part in re.split(r"(\d+)", label)
        if part
    )


def count_supported_sv_records(vcf_path):
    """Count supported SV records per chromosome; each BND line counts once."""
    counts = defaultdict(Counter)

    with pysam.VariantFile(os.fspath(vcf_path)) as vcf:
        for record in vcf:
            svtype = record.info.get("SVTYPE")
            if svtype in SV_TYPES:
                counts[record.chrom][svtype] += 1

    if not counts:
        raise ValueError(
            "The VCF contains no INS, DEL, INV, DUP, or BND records to plot."
        )

    return {
        chromosome: counts[chromosome]
        for chromosome in sorted(counts, key=natural_chromosome_key)
    }


def filter_small_contigs(
    counts, vcf_path, min_contig_fraction=0.001, visualize_small_contigs=False
):
    """Filter plotted contigs using a fraction of the total VCF genome length."""
    if visualize_small_contigs:
        return counts, "Including all contigs (--visualize-small-contigs)."

    with pysam.VariantFile(os.fspath(vcf_path)) as vcf:
        contig_lengths = {
            name: contig.length for name, contig in vcf.header.contigs.items()
        }

    # A partial genome length would produce an unreliable relative threshold.
    if not contig_lengths or any(length is None for length in contig_lengths.values()):
        return counts, (
            "Small-contig filtering was not applied because the VCF header does "
            "not provide lengths for every contig."
        )

    total_genome_length = sum(contig_lengths.values())
    minimum_length = total_genome_length * min_contig_fraction
    filtered_counts = {
        chromosome: chromosome_counts
        for chromosome, chromosome_counts in counts.items()
        if contig_lengths.get(chromosome, minimum_length) >= minimum_length
    }
    excluded_count = len(counts) - len(filtered_counts)

    if not filtered_counts:
        raise ValueError(
            "All contigs containing supported SV records fall below the computed "
            "length threshold. Use --visualize-small-contigs to include them."
        )

    message = (
        f"Genome length: {total_genome_length:,} bp. Excluded {excluded_count} "
        f"contig(s) shorter than {minimum_length:,.0f} bp "
        f"({min_contig_fraction * 100:g}% of genome length). "
        "Use --visualize-small-contigs to include them."
    )
    return filtered_counts, message


def plot_sv_counts(counts, output_path, input_name):
    """Write a stacked bar plot of SV record counts per chromosome."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib. Install it with 'pip install matplotlib'."
        ) from exc

    chromosomes = list(counts)
    figure_width = max(8, len(chromosomes) * 0.65)
    figure, axis = plt.subplots(figsize=(figure_width, 6))
    bottoms = [0] * len(chromosomes)

    for svtype in SV_TYPES:
        values = [counts[chromosome].get(svtype, 0) for chromosome in chromosomes]
        if not any(values):
            continue

        axis.bar(
            chromosomes,
            values,
            bottom=bottoms,
            color=SV_COLORS[svtype],
            label=svtype,
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]

    axis.set_title(f"Structural variant records per chromosome: {input_name}")
    axis.set_xlabel("Chromosome")
    axis.set_ylabel("Number of VCF records")
    axis.yaxis.set_major_locator(MaxNLocator(integer=True))
    axis.legend(
        title="SV type", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1)
    )
    axis.grid(axis="y", linestyle="--", alpha=0.3)
    axis.set_axisbelow(True)
    plt.setp(axis.get_xticklabels(), rotation=45, ha="right")
    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create a stacked bar plot of INS, DEL, INV, DUP, and individual "
            "BND record counts per chromosome."
        )
    )
    parser.add_argument("vcf", type=Path, help="Input VCF or compressed VCF file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("sv_counts_by_chromosome.png"),
        help="Output plot path (default: sv_counts_by_chromosome.png).",
    )
    parser.add_argument(
        "--min-contig-fraction",
        type=float,
        default=0.001,
        help=(
            "Minimum contig length as a fraction of total genome length "
            "(default: 0.001)."
        ),
    )
    parser.add_argument(
        "--visualize-small-contigs",
        action="store_true",
        help="Include all contigs regardless of their length.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        counts = count_supported_sv_records(args.vcf)
        if not 0 <= args.min_contig_fraction < 1:
            raise ValueError("--min-contig-fraction must be between 0 and 1.")
        counts, filter_message = filter_small_contigs(
            counts,
            args.vcf,
            args.min_contig_fraction,
            args.visualize_small_contigs,
        )
        print(filter_message)
        plot_sv_counts(counts, args.output, args.vcf.name)
    except (OSError, ValueError, RuntimeError) as exc:
        raise SystemExit(f"Error: {exc}") from exc

    print(f"Plot written to {args.output}")


if __name__ == "__main__":
    main()
