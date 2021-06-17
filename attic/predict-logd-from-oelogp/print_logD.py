import random

import click
from openeye.oechem import OEMol, OEParseSmiles
from openeye.oemolprop import OEGetXLogP


@click.command()
@click.option(
    "--fuzz",
    is_flag=True,
    default=False,
    show_default=True,
    help="Random change logP value by +/- 10%",
)
@click.argument(
	"solute", 
	default="CCCCCCCCO",
)
@click.argument(
	"solventa", 
	default="O",
)
@click.argument(
	"solventb", 
	default="CCCCCCCCO",
)
def get_LogD(solute, solventa, solventb, fuzz):
    mol = OEMol()
    OEParseSmiles(mol, solute)
    logP = OEGetXLogP(mol)
    if fuzz:
        # randomly change logP value by +/- 10%
        fuzzed_logP = logP + random.uniform(-0.1, 0.1) * logP
        click.echo(fuzzed_logP)
    else:
        click.echo(logP)


if __name__ == "__main__":
    get_LogD()
