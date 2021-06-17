import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles


@click.command()
@click.option("--smiles", help="SMILES string")
@click.argument("smiles_arg", required=False, default=None)
def calc_mol_wt(smiles, smiles_arg):
    if not smiles:
        smiles = smiles_arg
    if not smiles:
        click.echo("Must pass SMILES either with --smiles or directly", err=True)
    print(ExactMolWt(MolFromSmiles(smiles)))
