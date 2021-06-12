from unittest.mock import Mock, patch

import dask.distributed as dd
import pytest
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

    assert fake_run_container.call_count == 3

    fake_run_container.assert_called_with(
        "docker.io/mmh42/calc-subtract:0.1", '[{"diff": 3.0}, {"diff": 3.0}]'
    )
