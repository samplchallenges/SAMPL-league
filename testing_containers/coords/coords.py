import os.path
import argparse

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem

CONFORMATION_KEY = "conformation"
MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"


def calc_coords(output_dir, molfile, smiles,):

    if not output_dir:
        output_dir = ""

    def _printinfo(mol):
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)
        OUTPUT_FILENAME = "output.mol"
        output_path = os.path.join(output_dir, OUTPUT_FILENAME)
        with open(output_path, "w") as fp:
            print(Chem.MolToMolBlock(mol2), file=fp)
        print(CONFORMATION_KEY, output_path)
        print(MOLW_KEY, ExactMolWt(mol2))
        print(ATOMCOUNT_KEY, mol2.GetNumAtoms())
        print(BONDCOUNT_KEY, mol2.GetNumBonds())

    if molfile:
        mol = Chem.MolFromMolFile(molfile)
        _printinfo(mol)
        return 0
    if not smiles:
        click.echo("Must pass SMILES with --smiles", err=True)
        return 1
    mol = Chem.MolFromSmiles(smiles)
    _printinfo(mol)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', help="Output Directory")
    parser.add_argument('--molfile', help="MOL File")
    parser.add_argument('--smiles', help="SMILES string")

    args = parser.parse_args()

    calc_coords(args.output_dir, args.molfile, args.smiles)



