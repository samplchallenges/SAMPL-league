import pytest
from django.core.exceptions import ValidationError


def test_input_value_clean_challenge_match(smiles_molw_config, benzene_from_mol):
    input_value = benzene_from_mol.inputvalue_set.first()

    input_value.value_type.challenge = smiles_molw_config.challenge
    with pytest.raises(ValidationError) as exc_info:
        input_value.clean()

    assert list(exc_info.value.error_dict.keys()) == ["input_element"]


def test_input_value_contenttype_match(smiles_type, benzene_from_mol):
    input_value = benzene_from_mol.inputvalue_set.first()

    input_value.content_type = smiles_type.content_type

    with pytest.raises(ValidationError) as exc_info:
        input_value.clean()

    assert any("Value object's content type" in msg for msg in exc_info.value.messages)
