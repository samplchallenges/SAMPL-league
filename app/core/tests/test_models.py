# pylint: disable=unused-argument, unused-variable
import os.path
import re
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from django.conf import settings
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ValidationError
from django.db import models as django_models

from core import models
from core.tests import mocktime

TEST_DATA_PATH = Path(__file__).parent / "data"


def test_value_registration():
    with pytest.raises(ValidationError):

        @models.register_value_model
        class InvalidValueModel(models.GenericValue):
            class Meta:
                app_label = "core.apps.CoreConfig"

    @models.register_value_model
    class CharValueModel(models.GenericValue):
        value = django_models.CharField(max_length=100)

        class Meta:
            app_label = "core.apps.CoreConfig"


def test_container(scoring_container):
    assert scoring_container.uri == "ghcr.io/robbason/score-coords:latest"

    scoring_container.tag = None

    assert scoring_container.uri == "ghcr.io/robbason/score-coords"



def test_file_value(input_elements, molfile_type):
    elem = input_elements[0]
    challenge = elem.challenge
    hellofile = "hello.txt"

    with tempfile.TemporaryDirectory() as tmpdir:
        hellopath = os.path.join(tmpdir, hellofile)
        with open(hellopath, "w", encoding="utf-8") as fp:
            fp.write("Hello world")
            fp.flush()
        file_value = models.FileValue.from_string(hellopath, challenge=challenge)
        file_value.save()
        expected_dirpath = Path(f"file_uploads/challenges/{challenge.id}")
        dirpath = Path(file_value.value.name).parent
        assert dirpath == expected_dirpath
        assert file_value.challenge == challenge


def test_input_element(input_elements, benzene_from_mol):
    elem = input_elements[0]
    kwargs, file_kwargs = elem.all_values()
    assert kwargs == {"smiles": "c1ccccc1"}
    assert file_kwargs == {}

    kwargs, file_kwargs = benzene_from_mol.all_values()
    assert kwargs == {}
    molfile = file_kwargs["molfile"]
    dirpath = Path(molfile).parent
    relative_path = dirpath.relative_to(settings.MEDIA_ROOT)
    expected_path = f"file_uploads/challenges/{benzene_from_mol.challenge_id}"
    assert relative_path == Path(expected_path)


@pytest.mark.parametrize(['mocknow', 'compareval'], [(mocktime.inactive_before,False),
                                                     (mocktime.active,True),
                                                     (mocktime.inactive_after,False)])
def test_challenge(challenge, mocknow, compareval):
    with patch("django.utils.timezone.now", mocknow):
        assert challenge.is_active() == compareval


def test_load_prediction_file(
    container_factory, submission_factory, submission_run_factory, benzene_from_mol
):
    challenge = benzene_from_mol.challenge
    coordsfile_type = models.ValueType(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="conformation",
        description="3D output MOL file",
    )
    container = container_factory(challenge, "robbason/score-coords", "latest")
    submission = submission_factory(container)
    submission_run = submission_run_factory(submission)
    evaluation = models.Evaluation.objects.create(
        submission_run=submission_run, input_element=benzene_from_mol
    )
    filename = "Conformer3D_CID_241.mdl"
    output_path = TEST_DATA_PATH / filename

    mock_field_file_save = Mock()

    def mock_save_side_effect(name, content, **kwargs):
        assert name == filename

    mock_field_file_save.side_effect = mock_save_side_effect

    with patch("django.db.models.fields.files.FieldFile.save", mock_field_file_save):

        models.Prediction.load_output(
            challenge, evaluation, coordsfile_type, output_path
        )
        assert mock_field_file_save.call_count == 1


def test_container_arg(draft_submission, custom_string_arg, custom_file_arg):
    container = draft_submission.container
    assert container.custom_args() == {"stringarg": "hello world"}
    filepath = os.path.join(settings.MEDIA_ROOT, "container_args/{}/{}/filearg/example.*.txt".format(container.user_id, container.id))
    assert re.match(filepath, container.custom_file_args()["filearg"])
