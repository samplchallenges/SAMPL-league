import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles, MolFromMolFile


@click.command()
@click.option("--molfile", help="MOL File")
@click.option("--smiles", help="SMILES string")
@click.argument("smiles_arg", required=False, default=None)
def calc_mol_wt(molfile, smiles, smiles_arg):
    if molfile:
        print(ExactMolWt(MolFromMolFile(molfile)))
        return 0
    if not smiles:
        smiles = smiles_arg
    if not smiles:
        click.echo("Must pass SMILES either with --smiles or directly", err=True)
        return 1
    print(ExactMolWt(MolFromSmiles(smiles)))
    return 0
