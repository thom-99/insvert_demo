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
