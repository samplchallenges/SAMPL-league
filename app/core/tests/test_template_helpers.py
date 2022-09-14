import os.path
import re

from core import template_helpers


def test_template_table(smiles_docking_config_and_func):
    smiles_docking_config, add_element = smiles_docking_config_and_func
    add_element("benzene", "c1ccccc1", 99.0)
    add_element("methane", "C", 17.0)

    element_table = template_helpers.ElementTable(smiles_docking_config.challenge)

    type_keys = [t.key for t in element_table.types]
    assert type_keys == ["protein_pdb", "smiles"]

    assert len(element_table.rows) == 2

    row = element_table.rows[0]
    assert row[0] == "benzene"
    assert re.match("5qcr.+pdb", os.path.basename(row[1].value.name))

    assert row[2].value == "c1ccccc1"

    assert re.match("--smiles c1ccccc1 --protein_pdb 5qcr.+pdb", row[3])
