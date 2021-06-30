import os.path
import tempfile
import pytest
from django.core.exceptions import ValidationError
from django.db import models as django_models

from core import models


def test_value_registration():
    with pytest.raises(ValidationError):

        @register_value_model
        class InvalidValueModel(GenericValue):
            class Meta:
                app_label = "core.apps.CoreConfig"

    @register_value_model
    class CharValueModel(GenericValue):
        value = django_models.CharField(max_length=100)

        class Meta:
            app_label = "core.apps.CoreConfig"


def test_container(scoring_container):
    assert scoring_container.uri == "docker.io/mmh42/calc-subtract:0.1"

    scoring_container.tag = None

    assert scoring_container.uri == "docker.io/mmh42/calc-subtract"


def test_file_value(input_elements, challenge, molfile_type):
    elem = input_elements[0]
    hellofile = "hello.txt"
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, hellofile), "w") as fp:
            fp.write("Hello world")
            fp.flush()
        file_value = models.FileValue.from_string(
            hellofile, output_dir=tmpdir, challenge=challenge
        )
        file_value.save()
        assert (
            file_value.value.name == f"file_uploads/challenges/{challenge.id}/hello.txt"
        )
        assert file_value.challenge == challenge
