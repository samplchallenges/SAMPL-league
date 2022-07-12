# pylint: disable=unused-argument, unused-variable

import re
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from django.conf import settings
from django.contrib.contenttypes.models import ContentType

from core import batching, models
from core.tests import mocktime

TEST_DATA_PATH = Path(__file__).parent.parent / "data"


def test_container(scoring_container):
    assert scoring_container.uri == "ghcr.io/megosato/score-coords:latest"

    scoring_container.tag = None

    assert scoring_container.uri == "ghcr.io/megosato/score-coords"


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


@pytest.mark.parametrize("batch", (False, True))
def test_load_prediction_file(
    container_factory,
    submission_factory,
    submission_run_factory,
    benzene_from_mol,
    batch,
):
    challenge = benzene_from_mol.challenge
    coordsfile_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="conformation",
        description="3D output MOL file",
        batch_method="fake",
    )
    container = container_factory(challenge, "megosato/score-coords", "latest")
    submission = submission_factory(container)
    submission_run = submission_run_factory(submission)
    evaluation = models.Evaluation.objects.create(
        submission_run=submission_run, input_element=benzene_from_mol
    )
    filename = "Conformer3D_CID_241.mdl"
    output_path = TEST_DATA_PATH / filename
    if batch:
        challenge.max_batch_size = 2
        challenge.save()
        batch_group = challenge.current_batch_group()
        batch = batch_group.inputbatch_set.first()
        batch_evaluation = batch.batchevaluation_set.create(
            submission_run=submission_run
        )
        batch_path = str(TEST_DATA_PATH / "batched_mols.sdf")

        def mock_batch_save_side_effect(name, content, **kwargs):
            from rdkit import Chem

            assert Chem.MolFromMolFile(content.name) is not None

        mock_batch_file_save = Mock()
        mock_batch_file_save.side_effect = mock_batch_save_side_effect

        class FakeBatcher:
            @classmethod
            def invert(cls, batchfile):
                for _, value in batching.SDFBatcher.invert(batchfile):
                    yield benzene_from_mol.id, value
                    break

        batching.BATCHERS["fake"] = FakeBatcher

        with patch(
            "django.db.models.fields.files.FieldFile.save", mock_batch_file_save
        ):
            models.Prediction.load_batch_output(
                challenge, batch_evaluation, coordsfile_type, batch_path
            )
            assert mock_batch_file_save.call_count == 1

    else:
        mock_field_file_save = Mock()

        def mock_save_side_effect(name, content, **kwargs):
            assert name == filename

        mock_field_file_save.side_effect = mock_save_side_effect

        with patch(
            "django.db.models.fields.files.FieldFile.save", mock_field_file_save
        ):

            models.Prediction.load_evaluation_output(
                challenge, evaluation, coordsfile_type, output_path
            )
            assert mock_field_file_save.call_count == 1


def test_container_arg(draft_submission, custom_string_arg, custom_file_arg):
    container = draft_submission.container
    assert container.custom_args() == {"stringarg": "hello world"}
    filepath = str(
        Path(settings.MEDIA_ROOT)
        / "container_args/{}/{}/filearg/example.*.txt".format(
            container.user_id, container.id
        )
    )
    assert re.match(filepath, container.custom_file_args()["filearg"])


def test_cancel_requested(draft_submission, submission_run_factory):
    submission_run = submission_run_factory(draft_submission)
    submission_run.status = models.Status.CANCEL_PENDING
    submission_run.save()

    assert submission_run.check_cancel_requested()
