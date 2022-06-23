import csv
import os.path

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem #import MolFromSmiles, MolFromMolFile

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
        mol = Chem.MolFromMolFile(molfile)
        _printinfo(mol)
        return 0
    if not smiles:
        print("Must pass SMILES either with --smiles or directly")
        return 1
    mol = Chem.MolFromSmiles(smiles)
    _printinfo(mol)
    return 0


def batch_calc_mol_wt(output_dir, molfile, smiles):
    # It would be more efficient if all CSV outputs could share a file
    molw_file = os.path.join(output_dir, "molw.csv")
    atomc_file = os.path.join(output_dir, "atomc.csv")
    bondc_file = os.path.join(output_dir, "bondc.csv")
    molw_fp = open(molw_file, "w+")
    atomc_fp = open(atomc_file, "w+")
    bondc_fp = open(bondc_file, "w+")
    fieldnames = ["id", "name", "value"]
    molw_writer = csv.DictWriter(molw_fp, fieldnames=fieldnames)
    molw_writer.writeheader()

    atomc_writer = csv.DictWriter(atomc_fp, fieldnames=fieldnames)
    atomc_writer.writeheader()

    bondc_writer = csv.DictWriter(bondc_fp, fieldnames=fieldnames)
    bondc_writer.writeheader()


    def _printinfo(mol, id, name):
        row = {"id": id, "name": name}
        row["value"] = ExactMolWt(mol)
        molw_writer.writerow(row)
        row["value"] = mol.GetNumAtoms()
        atomc_writer.writerow(row)
        row["value"] = mol.GetNumBonds()
        bondc_writer.writerow(row)

    if molfile:
        # TODO: why did this fail as a context mgr?
        mol_supplier = Chem.SDMolSupplier(molfile)
        for idx, mol in enumerate(mol_supplier):
            _printinfo(mol, mol.GetIntProp("SAMPL_ID"), mol.GetProp("SAMPL_NAME"))
        print("Loaded", (idx+1), "mols from", molfile)
    else:
        if not smiles:
            print("Must pass SMILES either with --smiles or directly")
            return 1
        with open(smiles, "r") as smiles_fp:
            smiles_reader = csv.DictReader(smiles_fp)
            for row in smiles_reader:
                mol = Chem.MolFromSmiles(row["smiles"])
                _printinfo(mol, row["id"], row["name"])
    print(MOLW_KEY, "molw.csv")
    print(ATOMCOUNT_KEY, "atomc.csv")
    print(BONDCOUNT_KEY, "bondc.csv")
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--molfile', help="MOL File")
    parser.add_argument('--smiles', help="SMILES string")
    parser.add_argument('--output-dir', help="Output Dir")
    parser.add_argument('--batch', action="store_true", default=False, help="Input and output will be batch format")
    args = parser.parse_args()
    if args.batch:
        print("output dir", args.output_dir, os.path.exists(args.output_dir))
        print("molfile", args.molfile, os.path.exists(args.molfile))
        with open(args.molfile, "r") as fp:
            print(fp.read())
        batch_calc_mol_wt(args.output_dir, args.molfile, args.smiles)
    else:
        calc_mol_wt(args.molfile, args.smiles)
