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
from django.db import IntegrityError
from django.db import models as django_models

from core import models
from core.models import values
from core.tests import mocktime

TEST_DATA_PATH = Path(__file__).parent / "data"

# pylint: disable=protected-access


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


def test_container(scoring_container):
    assert scoring_container.uri == "ghcr.io/megosato/score-coords:latest"

    scoring_container.tag = None

    assert scoring_container.uri == "ghcr.io/megosato/score-coords"


def test_file_value(input_elements, molfile_type):
    elem = input_elements[0]
    challenge = elem.challenge
    hellofile = "hello.txt"

    with tempfile.TemporaryDirectory() as tmpdir:
        hellopath = os.path.join(tmpdir, hellofile)
        with open(hellopath, "w", encoding="utf-8") as fp:
            fp.write("Hello world")
            fp.flush()
        file_value = models.FileValue.from_string(
            hellopath, challenge=challenge, input_element=elem
        )
        file_value.save()
        expected_dirpath = Path(f"file_uploads/challenges/{challenge.id}")
        dirpath = Path(file_value.value.name).parent
        assert dirpath == expected_dirpath
        assert file_value.challenge == challenge


@pytest.mark.parametrize(
    ["mocknow", "compareval"],
    [
        (mocktime.inactive_before, False),
        (mocktime.active, True),
        (mocktime.inactive_after, False),
    ],
)
def test_challenge(challenge, mocknow, compareval):
    with patch("django.utils.timezone.now", mocknow):
        assert challenge.is_active() == compareval


def test_load_prediction_file(
    container_factory, submission_factory, submission_run_factory, benzene_from_mol
):
    challenge = benzene_from_mol.challenge
    coordsfile_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="conformation",
        description="3D output MOL file",
    )
    container = container_factory(challenge, "megosato/score-coords", "latest")
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

        models.Prediction.load_evaluation_output(
            challenge, evaluation, coordsfile_type, output_path
        )
        assert mock_field_file_save.call_count == 1


def test_container_arg(draft_submission, custom_string_arg, custom_file_arg):
    container = draft_submission.container
    assert container.custom_args() == {"stringarg": "hello world"}
    filepath = os.path.join(
        settings.MEDIA_ROOT,
        "container_args/{}/{}/filearg/example.*.txt".format(
            container.user_id, container.id
        ),
    )
    assert re.match(filepath, container.custom_file_args()["filearg"])


def test_cancel_requested(draft_submission, submission_run_factory):
    submission_run = submission_run_factory(draft_submission)
    submission_run.status = models.Status.CANCEL_PENDING
    submission_run.save()

    assert submission_run.check_cancel_requested()


def test_batch_build_works(smiles_molw_config, input_elements):
    group1 = models.InputBatchGroup.objects.create(
        challenge=smiles_molw_config.challenge, max_batch_size=1
    )
    batch1_priv = models.InputBatch.objects.create(batch_group=group1, is_public=False)
    batch1_pub = models.InputBatch.objects.create(batch_group=group1, is_public=True)

    public_elements = [element for element in input_elements if element.is_public]
    private_elements = [element for element in input_elements if not element.is_public]

    batch1_pub.batchup(public_elements)
    batch1_priv.batchup(private_elements)


def test_batch_checks(smiles_molw_config, input_elements):
    element = input_elements[1]
    assert element.is_public
    group = models.InputBatchGroup.objects.create(
        challenge=smiles_molw_config.challenge, max_batch_size=1
    )
    batch_priv = models.InputBatch.objects.create(batch_group=group)
    batch1_pub = models.InputBatch.objects.create(batch_group=group, is_public=True)
    batch2_pub = models.InputBatch.objects.create(batch_group=group, is_public=True)
    with pytest.raises(ValidationError):
        batch_priv._set_elements([element])

    batch1_pub._set_elements([element])

    with pytest.raises(IntegrityError):
        batch2_pub._set_elements([element])
