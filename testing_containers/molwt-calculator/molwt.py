import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles, MolFromMolFile

MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"


@click.command()
@click.option("--molfile", help="MOL File")
@click.option("--smiles", help="SMILES string")
@click.argument("smiles_arg", required=False, default=None)
def calc_mol_wt(molfile, smiles, smiles_arg):

    def _printinfo(mol):
        print(MOLW_KEY, ExactMolWt(mol))
        print(ATOMCOUNT_KEY, mol.GetNumAtoms())
        print(BONDCOUNT_KEY, mol.GetNumBonds())

    if molfile:
        mol = MolFromMolFile(molfile)
        _printinfo(mol)
        return 0
    if not smiles:
        smiles = smiles_arg
    if not smiles:
        click.echo("Must pass SMILES either with --smiles or directly", err=True)
        return 1
    mol = MolFromSmiles(smiles)
    _printinfo(mol)
    return 0
