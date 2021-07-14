import random

import click
from openeye.oechem import OEMol, OEParseSmiles
from openeye.oemolprop import OEGetXLogP
from rdkit import Chem
from rdkit.Chem import Crippen

@click.command()
@click.option(
    "--fuzz",
    is_flag=True,
    default=False,
    show_default=True,
    help="Random change logP value by +/- 10%",
)
@click.option(
    "--solute",
    default="CCCCCCCCO",
    help="solute SMILES string"
)
@click.option(
    "--solventa",
    default="O",
    help="solventA SMILES string"
)
@click.option(
    "--solventb",
    default="CCCCCCCCO",
    help="solventB SMILES string"
)
@click.option(
    "--output-dir", 
    help="Output Directory", 
    type=click.Path(exists=True)
)

def get_logd(solute, solventa, solventb, fuzz, output_dir):
    ''' takes in all inputs required for a LogD calculation (solute, solventA and
        solventB) but only calculates the LogP using the average of oechem and rdkit 
        logP. Ignores the solventA and solventB inputs
    '''
    oemol = OEMol()
    OEParseSmiles(oemol, solute)
    oelogD = OEGetXLogP(oemol)


    rdmol = Chem.MolFromSmiles(solute)
    rdlogD = Crippen.MolLogP(rdmol)

    logD = (oelogD + rdlogD)/2
    if fuzz:
        # randomly change logP value by +/- 10%
        fuzzed_logD = logD + random.uniform(-0.1, 0.1) * logD
        logD = fuzzed_logD
    print(f"LogD {logD})", end="")


if __name__ == "__main__":
    get_logd()
