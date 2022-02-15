from unittest.mock import Mock, patch

import dask.distributed as dd
import pytest
from django.contrib.contenttypes.models import ContentType
from django.core.management import call_command
from django.db import transaction

from core import models
from referee import scoring, tasks


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


def _run_and_check_evaluation(submission_run, evaluation):

    delayed = tasks.run_evaluation(
        submission_run.submission.id,
        evaluation.id,
        submission_run.id,
        True,
    )
    delayed.compute(scheduler="synchronous")
    assert submission_run.evaluation_set.count() == 1
    evaluation = submission_run.evaluation_set.get()
    assert evaluation.status == models.Status.SUCCESS
    prediction = evaluation.prediction_set.get()
    return prediction


def test_run_element_mol(molfile_molw_config, benzene_from_mol):
    submission_run = molfile_molw_config.submission_run
    evaluation = models.Evaluation.objects.create(
        input_element=benzene_from_mol, submission_run=submission_run
    )

    prediction = _run_and_check_evaluation(submission_run, evaluation)
    assert prediction.value == pytest.approx(78.046950192)


def test_run_element_custom(
    molfile_molw_config, benzene_from_mol, container_arg_factory
):
    submission_run = molfile_molw_config.submission_run
    submission = submission_run.submission
    container = submission.container
    container_arg_factory(container, key="stringarg", string_value="hello world")
    container_arg_factory(
        container, key="filearg", file_name="example.txt", file_body="Some text"
    )
    evaluation = models.Evaluation.objects.create(
        input_element=benzene_from_mol, submission_run=submission_run
    )
    with pytest.raises(AssertionError) as excinfo:
        _run_and_check_evaluation(submission_run, evaluation)

    assert "FAILURE" in str(excinfo.value)
    evaluation.refresh_from_db()
    assert "Error: No such option" in evaluation.log_stderr


def test_evaluation_scoring_failure(molfile_molw_config, benzene_from_mol):
    submission_run = molfile_molw_config.submission_run
    evaluation = models.Evaluation.objects.create(
        input_element=benzene_from_mol, submission_run=submission_run
    )

    with patch("referee.scoring.score_evaluation") as mock_score:
        mock_score.side_effect = scoring.ScoringFailureException("Mock failure")
        delayed = tasks.run_evaluation(
            submission_run.submission.id,
            evaluation.id,
            submission_run.id,
            True,
        )
        delayed.compute(scheduler="synchronous")

    evaluation.refresh_from_db()
    assert evaluation.log_stderr.endswith("Mock failure")
    assert evaluation.status == models.Status.FAILURE


@pytest.fixture
def evaluation_scores(smiles_molw_config, evaluations):
    score_types = smiles_molw_config.challenge.score_types[
        models.ScoreType.Level.EVALUATION
    ]

    def _score(evaluation):
        for score_type in score_types.values():
            yield models.EvaluationScore.create(
                score_type=score_type, evaluation=evaluation, value=99.9
            )

    return [_score(evaluation) for evaluation in evaluations]


def test_submission_run_scoring_failure(
    smiles_molw_config, evaluations, evaluation_scores
):
    submission_run = smiles_molw_config.submission_run
    with patch("referee.scoring._score_submission_run") as mock_score:
        mock_score.side_effect = scoring.ScoringFailureException("Mock failure")
        with pytest.raises(scoring.ScoringFailureException):
            delayed = tasks.check_and_score(
                submission_run.id, True, [models.Status.SUCCESS] * len(evaluations)
            )
            delayed.compute(scheduler="synchronous")

    submission_run.refresh_from_db()
    assert submission_run.log_stderr == "Mock failure"


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


def test_run_files(
    file_container, elem_factory, file_answer_key_factory, float_answer_key_factory
):
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
    molweight_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FloatValue),
        key="molWeight",
        description="Molecular Weight",
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
    molweight_answer = float_answer_key_factory(
        challenge, benzene_from_mol, molweight_type, 78.046950192
    )

    evaluation = models.Evaluation.objects.create(
        input_element=benzene_from_mol, submission_run=submission_run
    )
    delayed = tasks.run_evaluation(
        submission_run.submission.id,
        evaluation.id,
        submission_run.id,
        True,
    )
    delayed.compute(scheduler="synchronous")
    assert submission_run.evaluation_set.count() == 1
    evaluation = submission_run.evaluation_set.get()

    print(submission_run.evaluation_set.all())
    assert evaluation.status == models.Status.SUCCESS
    # prediction = evaluation.prediction_set.get()
    print(evaluation.prediction_set.all())
    assert prediction.value == pytest.approx(78.046950192)


def test_cancel_evaluation_before_run(molfile_molw_config, benzene_from_mol):
    submission_run = molfile_molw_config.submission_run
    evaluation = models.Evaluation.objects.create(
        input_element=benzene_from_mol, submission_run=submission_run
    )
    submission_run.mark_for_cancel()
    delayed = tasks.run_evaluation(
        submission_run.submission.id,
        evaluation.id,
        submission_run.id,
        True,
    )
    result = delayed.compute(scheduler="synchronous")
    assert result == models.Status.CANCELLED
    evaluation.refresh_from_db()
    assert evaluation.status == models.Status.CANCELLED


def test_cancel_submission_before_run(molfile_molw_config, benzene_from_mol):
    submission = molfile_molw_config.submission_run.submission
    delayed_conditional = tasks._trigger_submission_run(
        submission, True, is_public=True
    )
    submission.last_public_run().mark_for_cancel()
    result = delayed_conditional.compute(scheduler="synchronous")
    assert result is False
    assert submission.last_public_run().status == models.Status.CANCELLED
