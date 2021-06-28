import os.path

import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem

CONFORMATION_KEY = "conformation"
MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"


@click.command()
@click.option("--output-dir", help="Output Directory")
@click.option("--molfile", help="MOL File")
@click.option("--smiles", help="SMILES string")
@click.argument("smiles_arg", required=False, default=None)
def calc_coords(output_dir, molfile, smiles, smiles_arg):

    if not output_dir:
        output_dir = ""

    def _printinfo(mol):
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)
        OUTPUT_FILENAME = "output.mol"
        with open(os.path.join(output_dir, OUTPUT_FILENAME), "w") as fp:
            print(Chem.MolToMolBlock(mol2),file=fp)
        print(CONFORMATION_KEY, OUTPUT_FILENAME)
        print(MOLW_KEY, ExactMolWt(mol2))
        print(ATOMCOUNT_KEY, mol2.GetNumAtoms())
        print(BONDCOUNT_KEY, mol2.GetNumBonds())

    if molfile:
        mol = Chem.MolFromMolFile(molfile)
        _printinfo(mol)
        return 0
    if not smiles:
        smiles = smiles_arg
    if not smiles:
        click.echo("Must pass SMILES either with --smiles or directly", err=True)
        return 1
    mol = Chem.MolFromSmiles(smiles)
    _printinfo(mol)
    return 0
