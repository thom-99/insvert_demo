#!/usr/bin/env python3
"""Plot VCF variant positions across length-proportional chromosome spans."""

import argparse
import os
import re
from collections import defaultdict
from pathlib import Path

import pysam


VARIANT_TYPES = ("INS", "DEL", "INV", "DUP", "BND")
VARIANT_COLORS = {
    "INS": "#F2CF5B",
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


def infer_variant_type(record):
    """Read a supported structural-variant type from the record."""
    variant_type = record.info.get("SVTYPE") or record.info.get("VT")
    if variant_type in VARIANT_TYPES:
        return variant_type
    return None


def load_variant_positions(
    vcf_path, min_contig_fraction=0.001, visualize_small_contigs=False
):
    """Load supported variant positions and apply adaptive contig filtering."""
    positions = defaultdict(lambda: defaultdict(list))

    with pysam.VariantFile(os.fspath(vcf_path)) as vcf:
        contig_lengths = {
            name: contig.length for name, contig in vcf.header.contigs.items()
        }
        for record in vcf:
            variant_type = infer_variant_type(record)
            if variant_type is not None:
                positions[record.chrom][variant_type].append(record.pos)

    if not positions:
        raise ValueError("The VCF contains no supported variant records to plot.")

    plotted_chromosomes = set(positions)
    if visualize_small_contigs:
        filter_message = "Including all contigs (--visualize-small-contigs)."
    elif (
        not contig_lengths
        or any(length is None for length in contig_lengths.values())
        or any(chromosome not in contig_lengths for chromosome in plotted_chromosomes)
    ):
        filter_message = (
            "Small-contig filtering was not applied because the VCF header does "
            "not provide lengths for every contig."
        )
    else:
        total_genome_length = sum(contig_lengths.values())
        minimum_length = total_genome_length * min_contig_fraction
        retained_chromosomes = {
            chromosome
            for chromosome in plotted_chromosomes
            if contig_lengths[chromosome] >= minimum_length
        }
        excluded_count = len(plotted_chromosomes) - len(retained_chromosomes)
        positions = defaultdict(
            lambda: defaultdict(list),
            {
                chromosome: positions[chromosome]
                for chromosome in retained_chromosomes
            },
        )

        if not positions:
            raise ValueError(
                "All contigs containing supported variants fall below the computed "
                "length threshold. Use --visualize-small-contigs to include them."
            )

        filter_message = (
            f"Genome length: {total_genome_length:,} bp. Excluded {excluded_count} "
            f"contig(s) shorter than {minimum_length:,.0f} bp "
            f"({min_contig_fraction * 100:g}% of genome length). "
            "Use --visualize-small-contigs to include them."
        )

    missing_lengths = [
        chromosome for chromosome in positions if contig_lengths.get(chromosome) is None
    ]
    if missing_lengths:
        # Without lengths, proportional chromosome spans cannot be calculated.
        raise ValueError(
            "Cannot scale chromosome spans because length information is missing for: "
            + ", ".join(sorted(missing_lengths, key=natural_chromosome_key))
        )

    ordered_positions = {
        chromosome: positions[chromosome]
        for chromosome in sorted(positions, key=natural_chromosome_key)
    }
    return ordered_positions, contig_lengths, filter_message


def plot_variant_positions(
    positions, contig_lengths, output_path, input_name, chromosomes_per_row=4
):
    """Write a wrapped, length-proportional variant position plot."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator, StrMethodFormatter
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib. Install it with 'pip install matplotlib'."
        ) from exc

    chromosomes = list(positions)
    chromosome_rows = [
        chromosomes[index:index + chromosomes_per_row]
        for index in range(0, len(chromosomes), chromosomes_per_row)
    ]
    figure = plt.figure(
        figsize=(14, max(4.5, len(chromosome_rows) * 3.4)),
        constrained_layout=True,
    )
    outer_grid = figure.add_gridspec(len(chromosome_rows), 1, hspace=0.35)

    for row_index, chromosome_row in enumerate(chromosome_rows):
        width_ratios = [contig_lengths[chromosome] for chromosome in chromosome_row]

        row_grid = outer_grid[row_index].subgridspec(
            1, len(width_ratios), width_ratios=width_ratios, wspace=0.08
        )

        for column_index, chromosome in enumerate(chromosome_row):
            axis = figure.add_subplot(row_grid[0, column_index])
            for type_row, variant_type in enumerate(VARIANT_TYPES):
                x_values = positions[chromosome].get(variant_type, [])
                if not x_values:
                    continue

                # Deterministic vertical jitter separates collisions while
                # leaving every genomic X position unchanged.
                y_values = [
                    type_row
                    + (
                        ((position * 2654435761 + index * 1013904223) & 0xFFFFFFFF)
                        / 0xFFFFFFFF
                        - 0.5
                    )
                    * 0.5
                    for index, position in enumerate(x_values)
                ]
                axis.scatter(
                    x_values,
                    y_values,
                    s=7,
                    alpha=0.8,
                    color=VARIANT_COLORS[variant_type],
                    edgecolors="none",
                    rasterized=True,
                )

            axis.set_xlim(0, contig_lengths[chromosome])
            axis.set_ylim(-0.5, len(VARIANT_TYPES) - 0.5)
            axis.set_title(chromosome, fontsize=10)
            axis.set_yticks(range(len(VARIANT_TYPES)))
            if column_index == 0:
                axis.set_yticklabels(VARIANT_TYPES)
                axis.set_ylabel("Variant type")
            else:
                axis.set_yticklabels([])
                axis.tick_params(axis="y", length=0)
            axis.grid(axis="y", linestyle="--", alpha=0.25)
            axis.set_axisbelow(True)
            axis.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))

            axis.xaxis.set_major_locator(MaxNLocator(nbins=2, integer=True))
    figure.suptitle(f"Variant positions across the genome: {input_name}")
    figure.supxlabel("Position within chromosome (bp); panel width proportional to length")
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot INS, DEL, INV, DUP, and BND positions across "
            "length-proportional chromosome spans."
        )
    )
    parser.add_argument("vcf", type=Path, help="Input VCF or compressed VCF file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("variant_positions_by_chromosome.png"),
        help="Output plot path (default: variant_positions_by_chromosome.png).",
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
    parser.add_argument(
        "--chromosomes-per-row",
        type=int,
        default=4,
        help="Maximum chromosome panels per row (default: 4).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        if not 0 <= args.min_contig_fraction < 1:
            raise ValueError("--min-contig-fraction must be between 0 and 1.")
        if args.chromosomes_per_row < 1:
            raise ValueError("--chromosomes-per-row must be at least 1.")
        positions, contig_lengths, filter_message = load_variant_positions(
            args.vcf,
            args.min_contig_fraction,
            args.visualize_small_contigs,
        )
        print(filter_message)
        plot_variant_positions(
            positions,
            contig_lengths,
            args.output,
            args.vcf.name,
            args.chromosomes_per_row,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        raise SystemExit(f"Error: {exc}") from exc

    print(f"Plot written to {args.output}")


if __name__ == "__main__":
    main()
