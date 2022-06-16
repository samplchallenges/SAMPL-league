"""
Batch up file types.
"""
import csv
import os.path
import tempfile
from abc import ABC, abstractmethod

from rdkit import Chem

from . import values_helper

BATCHERS = {}


class BatcherBase(ABC):
    @staticmethod
    @abstractmethod
    def call(elements, key, output_path):
        pass

    @staticmethod
    @abstractmethod
    def invert(batchfile):
        pass


def register_batcher(cls):
    BATCHERS[cls.name] = cls
    return cls


@register_batcher
class CSVBatcher(BatcherBase):
    name = "csv"
    suffix = "csv"

    @staticmethod
    def call(elements, key, output_path):
        fieldnames = ("id", "name", "value")
        with open(output_path, "w", encoding="utf8") as csvfp:
            writer = csv.DictWriter(csvfp, fieldnames=fieldnames)
            writer.writeheader()
            for element in elements:
                values, _ = values_helper.all_values(element)
                row = {
                    "id": element.id,
                    "name": element.name,
                    "value": values[key],
                }
                writer.writerow(row)

    @staticmethod
    def invert(batchfile):
        with open(batchfile, encoding="utf8") as csvfp:
            reader = csv.DictReader(csvfp)
            for row in reader:
                yield row["id"], row["value"]


@register_batcher
class MolBatcher(BatcherBase):
    name = "mol"
    suffix = "mol"

    @staticmethod
    def call(elements, key, output_path):
        with Chem.SDWriter(output_path) as molwriter:
            for element in elements:
                _, file_values = values_helper.all_values(element)
                mol_file = file_values[key]
                mol = Chem.MolFromMolFile(mol_file)
                if mol is None:
                    raise Exception(f"Mol file {mol_file} could not be read")
                mol.SetIntProp("SAMPL_ID", element.id)
                mol.SetProp("SAMPL_NAME", element.name)
                molwriter.write(mol)

    @staticmethod
    def invert(batchfile):
        with tempfile.TemporaryDirectory() as temp_dir:
            with Chem.SDMolSupplier(batchfile) as mol_supplier:
                for mol in mol_supplier:
                    element_id = mol.GetIntProp("SAMPL_ID")
                    molpath = os.path.join(temp_dir, f"mol_{element_id}.sdf")
                    print(
                        Chem.MolToMolBlock(mol),
                        file=open(molpath, "w+", encoding="utf8"),
                    )
                    yield element_id, molpath
