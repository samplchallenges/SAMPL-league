import csv
import os.path
import tempfile

import pytest
from rdkit import Chem

from .. import batching, models

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_batch_csv(smiles_molw_config, input_elements):
    with tempfile.TemporaryDirectory() as dirname:
        output_path = os.path.join(dirname, "smiles_list.csv")
        batching.BATCHERS["csv"].call(input_elements, "smiles", output_path)
        with open(output_path, encoding="utf8") as fp:
            reader = csv.DictReader(fp)
            by_name = {}
            for row in reader:
                by_name[row["name"]] = row["value"]
        for elem in input_elements:
            assert elem.name in by_name
            found_smiles = by_name[elem.name]
            elem_smiles = (
                elem.inputvalue_set.filter(value_type=smiles_molw_config.input_type)
                .first()
                .value
            )
            assert found_smiles == elem_smiles


def test_invert_csv():
    batch_file = os.path.join(TEST_DATA_DIR, "batched_smiles.csv")
    rows = [(id, smiles) for id, smiles in batching.BATCHERS["csv"].invert(batch_file)]
    assert len(rows) == 2
    assert rows[0][0] == 1


def test_batch_mol(benzene_from_mol):
    with tempfile.TemporaryDirectory() as dirname:
        output_path = os.path.join(dirname, "mols.sdf")
        batching.BATCHERS["sdf"].call([benzene_from_mol], "molfile", output_path)

        suppl = Chem.SDMolSupplier(output_path)

        for idx, mol in enumerate(suppl):
            assert mol is not None
        # pylint: disable=undefined-loop-variable
        assert idx == 0
        assert mol.GetProp("SAMPL_ID") == str(benzene_from_mol.id)
        assert mol.GetProp("SAMPL_NAME") == benzene_from_mol.name


def test_invert_mol():
    batch_file = os.path.join(TEST_DATA_DIR, "batched_mols.sdf")
    rows = [(id, mol) for id, mol in batching.BATCHERS["sdf"].invert(batch_file)]
    assert len(rows) == 2
    assert rows[0][0] == 1
    molpath = rows[0][1]
    with Chem.SDMolSupplier(molpath) as reader:
        mol = next(reader)
    assert mol.GetProp("SAMPL_NAME") == "benzene"


@pytest.fixture
def corrupt_mol(molfile_molw_config, elem_factory):
    return elem_factory(
        molfile_molw_config.challenge,
        molfile_molw_config.input_type,
        "corrupt",
        "corrupt.mdl",
    )


def test_batch_invalid_mol(corrupt_mol):  # pylint: disable=redefined-outer-name
    with tempfile.TemporaryDirectory() as dirname:
        output_path = os.path.join(dirname, "mols.sdf")
        with pytest.raises(Exception):
            batching.BATCHERS["sdf"].call([corrupt_mol], "molfile", output_path)


def test_batch_hooks(
    molfile_molw_config, benzene_from_mol
):  # pylint: disable=unused-argument
    challenge = molfile_molw_config.challenge
    assert challenge.current_batch_group() is None

    challenge.max_batch_size = 10

    challenge.save()

    batch_group: models.InputBatchGroup = challenge.current_batch_group()

    assert batch_group is not None

    assert batch_group.inputbatch_set.count() == 1

    batch = batch_group.inputbatch_set.first()

    assert batch.batchfile_set.count() == 1
