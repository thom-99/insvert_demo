from .tests import input_files_tests
from . import simulate
from . import insert_streaming

import rich_click as click
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
click.rich_click.USE_RICH_MARKUP = True

# intializing Rich Console 
console = Console()

@click.group()
@click.version_option("0.0.0", prog_name="inSVert")
def cli():
    """
    [bold magenta]inSVert[/bold magenta]: Structural Variant Simulation & Insertion Toolkit.

    A tool to [bold green]simulate[/bold green] SVs into a VCF or [bold green]insert[/bold green] 
    existing variants into a FASTA reference genome.    

    • [bold green]simulate[/bold green] - Generate synthetic SVs and output them as VCF
    • [bold green]insert[/bold green] - Apply existing variants from VCF into a reference genome

    """
    pass


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
def simulate_cmd(config, reference, output):
    """
    Generate simulated structural variants based on a configuration file.
    
    This command creates synthetic SVs (insertions, deletions, duplications, 
    inversions) according to parameters specified in the CONFIG file and 
    outputs them in VCF format.
    
    [bold]Arguments:[/bold]
    
    CONFIG     - YAML configuration file specifying SV parameters

    REFERENCE  - Reference genome in FASTA format
    
    [bold]Example:[/bold]
    
      $ inSVert simulate my_config.yaml hg38.fasta -o my_variants.vcf
    """

    # 1. Header
    console.print(Panel(f"Running Simulation with [yellow]{config}[/yellow]", title="[bold cyan]inSVert Simulate[/bold cyan]", border_style="cyan"))

    # 2. Execution with Spinner
    with console.status("[bold cyan]Generating Structural Variants...[/bold cyan]", spinner="dots"):
        try:
            simulate.run(config, reference, output)
        except Exception as e:
            console.print(f"[bold red]Error:[/bold red] {e}")
            raise click.Abort()

    console.print(f"[bold green]✔ Done![/bold green] VCF written to [underline]{output}[/underline]\n")


# HANDLING THE INSERT MODULE 
@cli.command(name="insert")
@click.argument("reference", type=click.Path(exists=True, dir_okay=False))
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False))
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

def insert_cmd(reference, vcf, gc, output):
    """
    Insert structural variants from a VCF file into a reference genome.
    
    This command reads variants from a VCF file and applies them to the 
    reference genome, creating a modified FASTA file with the variants 
    incorporated.
    
    [bold]Arguments:[/bold]
    
    REFERENCE  - Reference genome in FASTA format
    
    VCF        - VCF file containing structural variants to insert
    
    [bold]Example:[/bold]
    
      $ inSVert insert hg38.fasta variants.vcf --gc 0.45 -o modified_genome.fasta
    """
    # 1. Header
    console.print(Panel(f"Inserting Variants from [yellow]{vcf}[/yellow]", title="[bold green]inSVert Insert[/bold green]", border_style="green"))    

    valid_vcf = input_files_tests.prepare_vcf(vcf)

    # 2. Execution with Spinner
    with console.status(f"[bold green]Processing Genome (GC={gc})...[/bold green]", spinner="dots"):
            try:
                insert_streaming.run(gc, reference, valid_vcf, output)
            except Exception as e:
                console.print(f"[bold red]Error:[/bold red] {e}")
                raise click.Abort()
                
    console.print(f"[bold green]✔ Done![/bold green] Modified genome saved to [underline]{output}[/underline]\n")
    
 






if __name__ == "__main__":
    cli()
