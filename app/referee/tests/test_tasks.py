
import dask.distributed as dd
import pytest
from django.contrib.contenttypes.models import ContentType
from django.core.management import call_command
from django.db import transaction

from core import models
from referee import tasks


@pytest.mark.django_db(transaction=True)
def test_run_and_score_submission():
    # This test will fail if run after another transaction=True test
    # See workaround in tests/test_views.py:test_run_submission

    transaction.commit()
    call_command("migrate", "core", "zero", interactive=False)
    call_command("migrate", "core", interactive=False)
    call_command("sample_data")
    transaction.commit()

    submission = models.Submission.objects.first()
    cluster = dd.LocalCluster(n_workers=4, preload=("daskworkerinit_tst.py",))
    dask_client = dd.Client(cluster)

    print(submission.id, submission)
    future = tasks.run_and_score_submission(dask_client, submission)
    result = future.result()
    assert result


def test_run_element_mol(molfile_molw_config, benzene_from_mol):
    submission_run = molfile_molw_config.submission_run
    delayed = tasks.run_element(
        submission_run.submission.id,
        benzene_from_mol.id,
        submission_run.id,
        True,
    )
    delayed.compute(scheduler="synchronous")
    assert submission_run.evaluation_set.count() == 1
    evaluation = submission_run.evaluation_set.get()
    assert evaluation.status == models.Status.SUCCESS
    prediction = evaluation.prediction_set.get()
    assert pytest.approx(prediction.value, 78.046950192)


@pytest.fixture
def file_container(challenge_factory, user, db):
    coord_challenge = challenge_factory("Coords Challenge")
    return models.Container.objects.create(
        name="CoordGenContainer",
        user=user,
        challenge=coord_challenge,
        registry="ghcr.io",
        label="robbason/calc-coords",
        tag="latest",
    )


def test_run_files(file_container, elem_factory, file_answer_key_factory):
    challenge = file_container.challenge
    scoring_container = models.Container.objects.create(
        name="subtraction container",
        user=file_container.user,
        challenge=challenge,
        registry="ghcr.io",
        label="robbason/score-coords",
        tag="latest",
    )
    score_maker = models.ScoreMaker.objects.create(
        challenge=challenge, container=scoring_container
    )
    models.ScoreType.objects.create(
        challenge=challenge, key="RMSE", level=models.ScoreType.Level.EVALUATION
    )
    models.ScoreType.objects.create(
        challenge=challenge, key="Avg RMSE", level=models.ScoreType.Level.SUBMISSION_RUN
    )
    molfile_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=True,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="molfile",
        description="2D input MOL file",
    )
    coordsfile_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="conformation",
        description="3D output MOL file",
    )

    submission = models.Submission.objects.create(
        name="Conformation Submission",
        user=file_container.user,
        container=file_container,
        challenge=challenge,
    )
    submission_run = models.SubmissionRun.objects.create(
        submission=submission,
        digest="cafef00d",
        is_public=True,
        status=models.Status.PENDING,
    )

    benzene_from_mol = elem_factory(
        challenge, molfile_type, "Benzene", "ChEBI_16716.mdl"
    )
    benzene_answer = file_answer_key_factory(
        challenge, benzene_from_mol, coordsfile_type, "Conformer3D_CID_241.mdl"
    )

    delayed = tasks.run_element(
        submission_run.submission.id,
        benzene_from_mol.id,
        submission_run.id,
        True,
    )
    delayed.compute(scheduler="synchronous")
    assert submission_run.evaluation_set.count() == 1
    evaluation = submission_run.evaluation_set.get()
    assert evaluation.status == models.Status.SUCCESS
    prediction = evaluation.prediction_set.get()
    assert pytest.approx(prediction.value, 78.046950192)
