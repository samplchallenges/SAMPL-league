# pylint: disable=unused-argument
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from django.contrib.contenttypes.models import ContentType

from core import batching, models

TEST_DATA_PATH = Path(__file__).parent.parent / "data"


@pytest.mark.parametrize("batch", (False, True))
def test_load_prediction_file(
    container_factory,
    submission_factory,
    submission_run_factory,
    file_answer_key_factory,
    benzene_from_mol,
    batch,
    caplog,
):
    import logging

    caplog.set_level(logging.DEBUG)
    challenge = benzene_from_mol.challenge
    coordsfile_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="conformation",
        description="3D output MOL file",
        batch_method="sdf",
    )
    container = container_factory(challenge, "megosato/score-coords", "latest")
    submission = submission_factory(container)
    submission_run = submission_run_factory(submission)
    evaluation = models.Evaluation.objects.create(
        submission_run=submission_run, input_element=benzene_from_mol
    )
    filename = "Conformer3D_CID_241.mdl"
    output_path = TEST_DATA_PATH / filename
    _akey = file_answer_key_factory(
        challenge, benzene_from_mol, coordsfile_type, output_path
    )
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
            # Need to explicitly delete the FileValue objects created here since
            # they will not have reasonable input_element_id values
            assert batch_evaluation.filevalue_set.count() == 2
            batch_evaluation.filevalue_set.all().delete()
            assert mock_batch_file_save.call_count == 2

        # Now verify can set values on parent
        pdbpose_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=False,
            on_parent_flag=True,
            content_type=ContentType.objects.get_for_model(models.FileValue),
            key="pdbpose",
            description="Prepared protein file PDB",
            batch_method="",
        )
        parent_elem = models.InputElement.objects.create(
            challenge=challenge,
            is_parent=True,
            name="Batch 1",
            is_public=True,
        )

        benzene_from_mol.parent = parent_elem
        benzene_from_mol.save()
        batch.parent_input_element = parent_elem
        batch.save()
        mock_batch_file_save.reset_mock(side_effect=None)
        with patch(
            "django.db.models.fields.files.FieldFile.save", mock_batch_file_save
        ):
            models.Prediction.load_batch_output(
                challenge, batch_evaluation, pdbpose_type, batch_path
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
