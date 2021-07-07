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
@click.option(
    "--smiles", 
    default="CCCCCCCCO"
)
@click.option(
    "--output-dir",
    type=click.Path(exists=True)
)
def get_LogP(smiles, fuzz, output_dir):
    mol = OEMol()
    OEParseSmiles(mol, smiles)
    logP = OEGetXLogP(mol)
    if fuzz:
        # randomly change logP value by +/- 10%
        fuzzed_logP = logP + random.uniform(-0.1, 0.1) * logP
        logP = fuzzed_logP
    print(f"LogP {logP}",end="")


if __name__ == "__main__":
    get_LogP()
