"""
Batch up file types.
"""
import csv
import os.path
import tempfile
from abc import ABC, abstractmethod

from rdkit import Chem

BATCHERS = {}


class BatcherBase(ABC):
    @classmethod
    @abstractmethod
    def call(cls, elements, key, output_path):
        pass


def register_batcher(cls):
    BATCHERS[cls.name] = cls
    return cls


@register_batcher
class CSVBatcher(BatcherBase):
    name = "csv"
    suffix = "csv"

    @classmethod
    def call(cls, elements, key, output_path):
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


@register_batcher
class MolBatcher(BatcherBase):
    name = "mol"
    suffix = "mol"

    @classmethod
    def call(cls, elements, key, output_path):
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


# TODO: transaction. In caller?
def generate_batches(challenge):
    group = challenge.inputbatchgroup_set.create(
        max_batch_size=challenge.max_batch_size
    )
    for is_public in (True, False):
        # TODO: each batch must share a parent as well
        elements = (
            challenge.inputelement_set.filter(is_public=is_public, is_parent=False)
            .order_by("name")
            .all()
        )
        for idx in range(0, len(elements), challenge.max_batch_size):
            batch_elements = elements[idx : idx + challenge.max_batch_size]
            batch = group.inputbatch_set.create(is_public=is_public)
            batchup(batch, batch_elements)


def batchup(batch, batch_elements):
    with tempfile.TemporaryDirectory() as temp_dir:
        for value_type in batch.batch_group.challenge.valuetype_set.filter(
            is_input_flag=True, on_parent_flag=False
        ):
            batcher = BATCHERS.get(value_type.batch_method)
            if batcher is None:
                raise ValueError(
                    f"Batch method must be set on value type {value_type.key}"
                )
            output_path = os.path.join(temp_dir, f"{value_type.key}.{batcher.suffix}")
            batcher.call(batch_elements, value_type.key, output_path)
    for batch_element in batch_elements:
        batch.add_element(batch_element)
