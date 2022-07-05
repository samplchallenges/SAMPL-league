import os.path

from core import template_helpers


def test_template_table(smiles_docking_config_and_func):
    smiles_docking_config, add_element = smiles_docking_config_and_func
    add_element("benzene", "c1ccccc1")
    add_element("methane", "C")

    element_table = template_helpers.ElementTable(smiles_docking_config.challenge)

    type_keys = [t.key for t in element_table.types]
    assert type_keys == ["protein_pdb", "smiles"]

    assert len(element_table.rows) == 2

    row = element_table.rows[0]
    assert row[0] == "benzene"
    assert os.path.basename(row[1].value.name) == "5qcr.pdb"

    assert row[2].value == "c1ccccc1"

    assert row[3] == "--smiles c1ccccc1 --protein_pdb 5qcr.pdb"
