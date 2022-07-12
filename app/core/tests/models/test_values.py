# pylint: disable=unused-argument, unused-variable

import os.path
import tempfile
from pathlib import Path

import pytest
from django.core.exceptions import ValidationError
from django.db import models as django_models

from core.models import values


def test_value_registration():
    with pytest.raises(ValidationError):

        @values.register_value_model
        class InvalidValueModel(values.GenericValue):
            class Meta:
                app_label = "core.apps.CoreConfig"

    @values.register_value_model
    class CharValueModel(values.GenericValue):
        value = django_models.CharField(max_length=100)

        class Meta:
            app_label = "core.apps.CoreConfig"


def test_file_value(input_elements, molfile_type):
    elem = input_elements[0]
    challenge = elem.challenge
    hellofile = "hello.txt"

    with tempfile.TemporaryDirectory() as tmpdir:
        hellopath = os.path.join(tmpdir, hellofile)
        with open(hellopath, "w", encoding="utf-8") as fp:
            fp.write("Hello world")
            fp.flush()
        file_value = values.FileValue.from_string(
            hellopath, challenge=challenge, input_element=elem
        )
        file_value.save()
        expected_dirpath = Path(f"file_uploads/challenges/{challenge.id}")
        dirpath = Path(file_value.value.name).parent
        assert dirpath == expected_dirpath
        assert file_value.challenge == challenge


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
