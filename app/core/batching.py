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
        "Assemble values for key into a file stored at output_path"

    @staticmethod
    @abstractmethod
    def invert(batchfile):
        "Given a file, yield tuples of ID, value"


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
                yield int(row["id"]), row["value"]


@register_batcher
class SDFBatcher(BatcherBase):
    name = "sdf"
    suffix = "sdf"

    @staticmethod
    def call(elements, key, output_path):
        with Chem.SDWriter(output_path) as writer:
            for element in elements:
                _, file_values = values_helper.all_values(element)
                sdfile = file_values[key]
                with Chem.SDMolSupplier(sdfile) as reader:
                    mol = next(reader)
                if mol is None:
                    raise Exception(f"SDF {sdfile} could not be read")
                mol.SetIntProp("SAMPL_ID", element.id)
                mol.SetProp("SAMPL_NAME", element.name)
                writer.write(mol)

    @staticmethod
    def invert(batchfile):
        temp_dir = tempfile.mkdtemp()
        with Chem.SDMolSupplier(batchfile) as mol_supplier:
            for mol in mol_supplier:
                element_id = mol.GetIntProp("SAMPL_ID")
                sdpath = os.path.join(temp_dir, f"mol_{element_id}.sdf")
                with Chem.SDWriter(sdpath) as writer:
                    writer.write(mol)
                yield element_id, sdpath
