"""
Batch up file types. Could be run from a container so we don't have to install RDKit into the main virtualenv (TBD)
"""
import csv

from rdkit import Chem

BATCHERS = {}


def _batcher(key):
    def func_register(func):
        BATCHERS[key] = func
        # TODO: verify function signature with inspect.getargspec
        return func

    return func_register


@_batcher("csv")
def batch_csv(elements, key, output_path):
    fieldnames = ("id", "name", "value")
    with open(output_path, "w", encoding="utf8") as csvfp:
        writer = csv.DictWriter(csvfp, fieldnames=fieldnames)
        writer.writeheader()
        for element in elements:
            values, _ = element.all_values()
            row = {
                "id": element.id,
                "name": element.name,
                "value": values[key],
            }
            writer.writerow(row)


@_batcher("mol")
def batch_mol(elements, key, output_path):
    with Chem.SDWriter(output_path) as molwriter:
        for element in elements:
            _, file_values = element.all_values()
            mol_file = file_values[key]
            mol = Chem.MolFromMolFile(mol_file)
            if mol is None:
                raise Exception(f"Mol file {mol_file} could not be read")
            mol.SetIntProp("SAMPL_ID", element.id)
            mol.SetProp("SAMPL_NAME", element.name)
            molwriter.write(mol)
