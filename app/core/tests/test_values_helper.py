from pathlib import Path

from django.conf import settings

from core import values_helper


def test_input_element(input_elements, benzene_from_mol):
    elem = input_elements[0]
    kwargs, file_kwargs = values_helper.all_values(elem)
    assert kwargs == {"smiles": "c1ccccc1"}
    assert file_kwargs == {}

    kwargs, file_kwargs = values_helper.all_values(benzene_from_mol)
    assert kwargs == {}
    molfile = file_kwargs["molfile"]
    dirpath = Path(molfile).parent
    relative_path = dirpath.relative_to(settings.MEDIA_ROOT)
    expected_path = f"file_uploads/challenges/{benzene_from_mol.challenge_id}"
    assert relative_path == Path(expected_path)


def test_parent_element(smiles_docking_config_and_func):
    _smiles_docking_config, add_element_func = smiles_docking_config_and_func
    benz = add_element_func("benzene", "c1ccccc1")
    _methane = add_element_func("methane", "C")

    _values, file_values = values_helper.all_values(benz)
    assert "protein_pdb" in file_values


def test_batch_values(smiles_docking_config_and_func):
    smiles_docking_config, add_element_func = smiles_docking_config_and_func
    challenge = smiles_docking_config.challenge
    mols = [("benzene", "c1ccccc1"), ("methane", "C"), ("water", "O")]
    for name, smiles in mols:
        add_element_func(name, smiles)
    challenge.max_batch_size = 3
    challenge.save()
    batch_group = challenge.current_batch_group()
    batch = batch_group.inputbatch_set.first()
    values, file_values = values_helper.batch_values(batch)
    assert values == {}
    assert file_values.keys() == {"protein_pdb", "smiles"}
