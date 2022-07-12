import os.path
import argparse
import csv
import sys

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem

CONFORMATION_KEY = "conformation"
MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"


def batch_coords(output_dir, sdfile, smilesfile):
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

    def _printinfo(mol, idx, name):
        row = {"id": idx, "name": name}
        row["value"] = ExactMolWt(mol)
        molw_writer.writerow(row)
        row["value"] = mol.GetNumAtoms()
        atomc_writer.writerow(row)
        row["value"] = mol.GetNumBonds()
        bondc_writer.writerow(row)

    output_file = os.path.join(output_dir, "coords.sdf")
    writer = Chem.SDWriter(output_file)

    for id, name, mol in _itermols(sdfile, smilesfile):
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)
        writer.write(mol2)
        _printinfo(mol, id, name)
    print(CONFORMATION_KEY, output_file)
    print(MOLW_KEY, "molw.csv")
    print(ATOMCOUNT_KEY, "atomc.csv")
    print(BONDCOUNT_KEY, "bondc.csv")


def _itermols(sdfile, smilesfile):
    if sdfile:
        mol_supplier = Chem.SDMolSupplier(sdfile)
        for mol in mol_supplier:
            yield mol.GetIntProp("SAMPL_ID"), mol.GetProp("SAMPL_NAME"), mol
    elif smilesfile:
        with open(smilesfile, "r") as smiles_fp:
            smiles_reader = csv.DictReader(smiles_fp)
            for row in smiles_reader:
                mol = Chem.MolFromSmiles(row["value"])
                yield row["id"], row["name"], mol
    else:
        raise ValueError("Either molfile or smilesfile must be set")


def calc_coords(
    output_dir,
    molfile,
    smiles,
):

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
        print("Must pass SMILES with --smiles", file=sys.stderr)
        return 1
    mol = Chem.MolFromSmiles(smiles)
    _printinfo(mol)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", help="Output Directory")
    parser.add_argument("--molfile", help="MOL File")
    parser.add_argument("--smiles", help="SMILES string")
    parser.add_argument("--protein_pdb", help="For testing")
    parser.add_argument(
        "--batch",
        action="store_true",
        default=False,
        help="Input and output will be batch format",
    )

    args = parser.parse_args()
    if args.protein_pdb:
        print("Protein PDB:", args.protein_pdb)
    if args.batch:
        batch_coords(args.output_dir, args.molfile, args.smiles)
    else:
        calc_coords(args.output_dir, args.molfile, args.smiles)
