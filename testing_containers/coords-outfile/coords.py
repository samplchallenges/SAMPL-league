import os.path
import argparse

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem

CONFORMATION_KEY = "conformation"
MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"
OUTFILE_KEY = "outfile"
def calc_coords(output_dir, molfile, smiles):

    if not output_dir:
        output_dir = ""

    def _printinfo(mol):
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)
        OUTPUT_FILENAME = "output.mol"
        output_path = os.path.join(output_dir, OUTPUT_FILENAME)
        with open(output_path, "w") as fp:
            print(Chem.MolToMolBlock(mol2), file=fp)
        print(MOLW_KEY, ExactMolWt(mol2))
        print(ATOMCOUNT_KEY, mol2.GetNumAtoms())
        print(BONDCOUNT_KEY, mol2.GetNumBonds())

    if molfile:
        mol = Chem.MolFromMolFile(molfile)
        _printinfo(mol)
        return 0
    if not smiles:
        print("Must pass SMILES either with --smiles or directly")
        return 1
    mol = Chem.MolFromSmiles(smiles)

    with open(os.path.join(output_dir, "outfile.txt"), 'w') as of:
        of.write(f"Keys: {MOLW_KEY}, {ATOMCOUNT_KEY}, {BONDCOUNT_KEY}") 
    _printinfo(mol)   
    print(OUTFILE_KEY, os.path.join(output_dir, "outfile.txt"))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--molfile")
    parser.add_argument("--smiles")
    parser.add_argument("--output-dir")

    args = parser.parse_args()
    calc_coords(args.output_dir, args.molfile, args.smiles)

