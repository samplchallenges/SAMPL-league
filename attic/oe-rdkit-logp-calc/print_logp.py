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
@click.argument("smiles", default="CCCCCCCCO")
def get_logp(smiles, fuzz):
    oemol = OEMol()
    OEParseSmiles(oemol, smiles)
    oelogP = OEGetXLogP(oemol)


    rdmol = Chem.MolFromSmiles(smiles)
    rdlogP = Crippen.MolLogP(rdmol)

    logP = (oelogP + rdlogP)/2

    print(logP)


if __name__ == "__main__":
    get_logp()
