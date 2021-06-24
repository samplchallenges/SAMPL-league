import click
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
	rdmol = Chem.MolFromSmiles(smiles)
	rdlogP = Crippen.MolLogP(rdmol)
	print(rdlogP)


if __name__ == "__main__":
	get_logp()
