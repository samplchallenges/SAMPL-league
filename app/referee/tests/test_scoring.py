# pylint: skip-file
import pytest

from core import models
from referee import scoring


def test_score_submission(
    smiles_molw_config,
    evaluations,
    evaluation_scores,
):
    submission_run = smiles_molw_config.submission_run
    submission = submission_run.submission

    scoring.score_submission(submission.pk, submission_run.pk)

    score_type = models.ScoreType.objects.get(
        challenge=smiles_molw_config.challenge, key="rmse"
    )

    assert models.SubmissionRunScore.objects.count() == 1
    submission_run_score = models.SubmissionRunScore.objects.get(score_type=score_type)
    assert submission_run_score.value == pytest.approx(3.0)

    # assert models.Prediction.objects.count() == 2
    assert models.EvaluationScore.objects.count() == 2

    evaluation = evaluations[0]
    assert evaluation.scores.count() == 1
    score = evaluation.scores.first()
    assert score.value == pytest.approx(3.0)


#    fake_run.assert_called_with(
#        "docker.io/mmh42/calc-subtract:0.1", '[{"diff": 3.0}, {"diff": 3.0}]'
#    )
