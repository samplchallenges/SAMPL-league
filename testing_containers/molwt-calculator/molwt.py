from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles, MolFromMolFile

import argparse

MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"


def calc_mol_wt(molfile, smiles):

    def _printinfo(mol):
        print(MOLW_KEY, ExactMolWt(mol))
        print(ATOMCOUNT_KEY, mol.GetNumAtoms())
        print(BONDCOUNT_KEY, mol.GetNumBonds())

    if molfile:
        mol = MolFromMolFile(molfile)
        _printinfo(mol)
        return 0
    if not smiles:
        click.echo("Must pass SMILES either with --smiles or directly", err=True)
        return 1
    mol = MolFromSmiles(smiles)
    _printinfo(mol)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--molfile', help="MOL File")
    parser.add_argument('--smiles', help="SMILES string")
    parser.add_argument('--output-dir', help="Output Dir")
    args = parser.parse_args()

    calc_mol_wt(args.molfile, args.smiles)
