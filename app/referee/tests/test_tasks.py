from unittest.mock import Mock, patch

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


def test_run_element_mol(submission_run_public, benzene_from_mol):
    delayed = tasks.run_element(
        submission_run_public.submission.id,
        benzene_from_mol.id,
        submission_run_public.id,
        True,
    )
    delayed.compute(scheduler="synchronous")
    assert submission_run_public.evaluation_set.count() == 1
    evaluation = submission_run_public.evaluation_set.get()
    assert evaluation.status == models.Status.SUCCESS
    prediction = evaluation.prediction_set.get()
    assert pytest.approx(prediction.value, 78.046950192)


@pytest.fixture
def submission_run_public(draft_submission, db):
    return models.SubmissionRun.objects.create(
        submission=draft_submission,
        digest="cafef00d",
        is_public=True,
        status=models.Status.PENDING,
    )


@pytest.fixture
def evaluations(challenge, submission_run_public, input_elements, molw_type, db):
    evaluations_list = []
    value = "97.08"
    for input_element in input_elements:
        if input_element.is_public:
            evaluation = models.Evaluation.objects.create(
                input_element=input_element, submission_run=submission_run_public
            )
            evaluations_list.append(evaluation)
            tasks._save_prediction(challenge, evaluation, molw_type, value)

    return evaluations_list


@pytest.fixture
def evaluation_scores(challenge, evaluations, score_types):
    score_value = 3.0
    evaluation_score_type = challenge.scoretype_set.filter(
        level=models.ScoreType.Level.EVALUATION
    ).get()

    return [
        models.EvaluationScore.objects.create(
            evaluation=evaluation, score_type=evaluation_score_type, value=score_value
        )
        for evaluation in evaluations
    ]


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


def test_save_prediction(challenge, submission_run_public, input_elements, molw_type):
    input_element = input_elements[0]
    evaluation = models.Evaluation.objects.create(
        input_element=input_element, submission_run=submission_run_public
    )
    value = "95.0"
    tasks._save_prediction(challenge, evaluation, molw_type, value)


def test_score_submission(
    draft_submission,
    submission_run_public,
    scoring_container,
    score_maker,
    evaluations,
    evaluation_scores,
    score_types,
):

    fake_run_container = Mock(return_value="3.0")
    mock_ever_given = Mock()
    mock_ever_given.run_container = fake_run_container
    with patch("ever_given.wrapper", mock_ever_given):
        tasks.score_submission(draft_submission.pk, submission_run_public.pk)

    assert models.SubmissionRunScore.objects.count() == 1
    submission_run_score = models.SubmissionRunScore.objects.get()
    assert submission_run_score.value == pytest.approx(3.0)

    assert models.Prediction.objects.count() == 2
    assert models.EvaluationScore.objects.count() == 2

    evaluation = evaluations[0]
    assert evaluation.scores.count() == 1
    score = evaluation.scores.first()
    assert score.value == pytest.approx(3.0)

    assert fake_run_container.call_count == 1

    fake_run_container.assert_called_with(
        "docker.io/mmh42/calc-subtract:0.1", '[{"diff": 3.0}, {"diff": 3.0}]'
    )
