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
    help="Randomly change logP value by +/- 10% [default: False]",
)
@click.option(
    "--solute", 
    default="CCCCCCCCO",
    help="solute SMILES string"
)
@click.option(
    "--solventA", 
    default="O",
    help="solventA SMILES string"
)
@click.option(
    "--solventB", 
    default="CCCCCCCCO",
    help="solventB SMILES string"
)
def get_logd(solute, solventA, solventB, fuzz):
    ''' takes in all inputs required for a LogD calculation (solute, solventA and
        solventB) but only calculates the LogP using oechem and ignores the solventA 
        and solventB inputs
    '''
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
    get_logd()
