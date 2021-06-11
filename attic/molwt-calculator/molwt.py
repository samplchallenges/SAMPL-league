import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles


@click.command()
@click.argument("smiles")
def calc_mol_wt(smiles):
    print(ExactMolWt(MolFromSmiles(smiles)))
