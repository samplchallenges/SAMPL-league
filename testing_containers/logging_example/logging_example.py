"""
To illustrate logging capabilities, runs slowly, printing out messages on
stdout and stderr
"""
import os.path
import sys
import time

import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem

CONFORMATION_KEY = "conformation"
MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"

def _p():
    time.sleep(1)

def mprint(*args, **kwargs):
    kwargs["flush"] = kwargs.get("flush", True)
    print(*args, **kwargs)


@click.command()
@click.option("--output-dir", help="Output Directory", type=click.Path(exists=True))
@click.option("--molfile", help="MOL File", type=click.Path(exists=True))
@click.option("--smiles", help="SMILES string")
@click.argument("smiles_arg", required=False, default=None)
def calc_coords(output_dir, molfile, smiles, smiles_arg):
    mprint("Starting coords")
    _p()
    mprint("Error message #1", file=sys.stderr)
    _p()
    for _x in range(5):
        mprint(".", end="")
        _p()
    mprint("Done initializing")
    mprint("Error message #2", file=sys.stderr)
    _p()
    if not output_dir:
        output_dir = ""

    def _printinfo(mol):
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)
        OUTPUT_FILENAME = "output.mol"
        output_path = os.path.join(output_dir, OUTPUT_FILENAME)
        with open(output_path, "w") as fp:
            mprint(Chem.MolToMolBlock(mol2), file=fp)
        mprint(CONFORMATION_KEY, output_path)
        _p()
        mprint(MOLW_KEY, ExactMolWt(mol2))
        _p()
        mprint(ATOMCOUNT_KEY, mol2.GetNumAtoms())
        _p()
        mprint(BONDCOUNT_KEY, mol2.GetNumBonds())
        _p()

    if molfile:
        mol = Chem.MolFromMolFile(molfile)
        _mprintinfo(mol)
        return 0
    if not smiles:
        smiles = smiles_arg
    if not smiles:
        click.echo("Must pass SMILES either with --smiles or directly", err=True)
        return 1
    mol = Chem.MolFromSmiles(smiles)
    mprint("Error message #3", file=sys.stderr)
    _printinfo(mol)
    return 0
