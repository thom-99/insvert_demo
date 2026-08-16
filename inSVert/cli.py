from . import input_validation
from . import simulate
from . import insert

import rich_click as click
from pathlib import Path
from textwrap import dedent
from rich.console import Console
from rich.panel import Panel
click.rich_click.USE_RICH_MARKUP = True

# intializing Rich Console 
console = Console()

@click.group()
@click.version_option("0.1.0", prog_name="inSVert")
def cli():
    """
    [bold magenta]inSVert[/bold magenta]: Structural Variant Simulation & Insertion Toolkit.

    A tool to [bold green]simulate[/bold green] SVs into a VCF or [bold green]insert[/bold green] 
    existing variants into a FASTA reference genome.    

    • [bold green]simulate[/bold green] - Generate synthetic SVs and output them as VCF
    • [bold green]insert[/bold green] - Apply existing variants from VCF into a reference genome

    """
    pass


@cli.command(name="generate-configfile")
@click.option(
    "--output",
    "-o",
    type=click.Path(dir_okay=False, path_type=Path),
    default=Path("config.yaml"),
    show_default=True,
    help="Path for the generated configuration file.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite the output file if it already exists.",
)
def generate_configfile(output, force):
    """Generate a template YAML configuration file."""

    if output.exists() and not force:
        raise click.ClickException(
            f"'{output}' already exists. Use --force to overwrite it."
        )

    template = dedent(
        """\
        # inSVert simulation configuration

        genome:
          ploidy: 2
          heterozygosity: 0.5

        variants:
          INS:
            count: 100
            distribution: normal
            parameters:
              median_length: 500
              min_length: 50
              max_length: 5000
              sigma: 100

          DEL:
            count: 100
            distribution: pareto
            parameters:
              median_length: 500
              min_length: 50
              max_length: 5000

          DUP:
            count: 50
            distribution: normal
            parameters:
              median_length: 1000
              min_length: 100
              max_length: 10000
            copy_number:
              min: 2
              max: 5
              weights: [0.5, 0.3, 0.15, 0.05]

          SNP:
            ratio: 0.0001
            tstv_ratio: 2.0
        """
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(template, encoding="utf-8")

    console.print(
        f"[bold green]✔ Configuration template written to "
        f"[underline]{output}[/underline][/bold green]"
    )


# HANDLING THE SIMULATE MODULE 
@cli.command(name="simulate")
@click.argument("config", type=click.Path(exists=True, dir_okay=False))
@click.argument("reference", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o", "--output", 
    default="simulated.vcf", 
    show_default=True,
    help="Path where the output VCF will be saved."
)
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Global random seed for reproducible simulations"
)
@click.option(
    "--exclude",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="BED file containing genomic coordinates to exclude from SV simulation."
)
@click.option(
    "--non-symbolic",
    is_flag=True,
    default=False,
    help="Output explicit REF/ALT sequences instead of symbolic tags (<INS>, <DEL>, etc.)."
)
def simulate_cmd(config, reference, output, seed, exclude, non_symbolic):
    """
    Generate simulated structural variants based on a configuration file.
    
    This command creates synthetic SVs (insertions, deletions, duplications, 
    inversions) according to parameters specified in the CONFIG file and 
    outputs them in VCF format.
    
    [bold]Arguments:[/bold]
    
    CONFIG     - YAML configuration file specifying SV parameters

    REFERENCE  - Reference genome in FASTA format
    
    [bold]Example:[/bold]
    
      $ inSVert simulate my_config.yaml hg38.fasta --seed 123 -o my_variants.vcf
    """
    # input validation
    console.print("[yellow]Validating input files...[/yellow]")
    try:
        input_validation.validate_fasta(reference)
        input_validation.validate_bed(exclude, reference)
        input_validation.validate_config(config)
        console.print("[bold green]✓ Inputs verified successfully![/bold green]")

    except Exception as e:
        console.print(f"\n[bold red]Validation Error:[/bold red] {e}")
        raise click.Abort
    
    # execution
    # 1. Header
    console.print(Panel(f"Running Simulation with [yellow]{config}[/yellow]", title="[bold cyan]inSVert Simulate[/bold cyan]", border_style="cyan"))

    # 2. Execution with Spinner
    with console.status("[bold cyan]Generating Structural Variants...[/bold cyan]", spinner="dots"):
        try:
            simulate.run(config, reference, output, seed, exclude, non_symbolic)
        except Exception as e:
            console.print(f"[bold red]Error:[/bold red] {e}")
            raise click.Abort()

    console.print(f"[bold green]✔ Done![/bold green] VCF written to [underline]{output}[/underline]\n")


# HANDLING THE INSERT MODULE 
@cli.command(name="insert")
@click.argument("reference", type=click.Path(exists=True, dir_okay=False))
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--ploidy", 
    type=click.IntRange(min=1), 
    required=True, 
    help="Ploidy of the simulated organism (number of haplotype copies to generate)."
)
@click.option(
    "--gc", 
    type=click.FloatRange(0.0, 1.0), 
    default=0.41, 
    show_default=True,
    help="Target GC content for random sequence generation."
)
@click.option(
    "-o", "--output", 
    default="edited_genome.fasta", 
    show_default=True, 
    help="Path for the modified output FASTA."
)
@click.option(
    "--skip-unphased-variants",
    is_flag=True,
    default=False,
    help="Skip variants with unphased genotypes (e.g., 0/1) instead of randomly assigning them to a haplotype."
)
@click.option(
    "--sample-name",
    default="Sample",
    show_default=True,
    help="Name used in multi-haplotype FASTA headers."
)
@click.option(
    "--split-haplotypes",
    is_flag=True,
    default=False,
    help="Write each haplotype to a separate FASTA file using the output name plus _hapN."
)

def insert_cmd(reference, vcf, ploidy, gc, output, skip_unphased_variants, sample_name, split_haplotypes):
    """
    Insert structural variants from a VCF file into a reference genome.
    
    This command reads variants from a VCF file and applies them to the 
    reference genome, creating a modified FASTA file with the variants 
    incorporated.
    
    [bold]Arguments:[/bold]
    
    REFERENCE  - Reference genome in FASTA format
    
    VCF        - VCF file containing structural variants to insert
    
    --ploidy   - the number of haplotype copies you want to produce

    [bold]Example:[/bold]
    
      $ inSVert insert hg38.fasta variants.vcf --ploidy 2 --gc 0.45 -o modified_genome.fasta
    """

    # input validation
    try:
        input_validation.validate_fasta(reference)
        input_validation.validate_vcf(vcf, reference, ploidy)
        input_validation.validate_output_path(reference,output)

    
    except Exception as e:
        console.print(f"\n[bold red]Validation Error:[/bold red] {e}")
        raise click.Abort()

    # execution
    # 1. Header
    console.print(Panel(f"Inserting Variants from [yellow]{vcf}[/yellow]", title="[bold green]inSVert Insert[/bold green]", border_style="green"))    

    # 2. Execution with Spinner
    with console.status(f"[bold green]Processing Genome (GC={gc})...[/bold green]", spinner="dots"):
            try:
                insert.run(gc, reference, vcf, ploidy, output, skip_unphased_variants, sample_name, split_haplotypes)
            except Exception as e:
                console.print(f"[bold red]Error:[/bold red] {e}")
                raise click.Abort()
                
    console.print(f"[bold green]✔ Done![/bold green] Modified genome saved to [underline]{output}[/underline]\n")
    
 






if __name__ == "__main__":
    cli()
