# pylint: skip-file
import pytest

from core import models


@pytest.fixture
def evaluations(smiles_molw_config, input_elements):
    challenge = smiles_molw_config.challenge
    submission_run = smiles_molw_config.submission_run
    molw_type = smiles_molw_config.output_type
    evaluations_list = []
    value = "97.08"
    for input_element in input_elements:
        if input_element.is_public:
            evaluation = models.Evaluation.objects.create(
                input_element=input_element, submission_run=submission_run
            )
            evaluations_list.append(evaluation)
            models.Prediction.load_evaluation_output(
                challenge, evaluation, molw_type, value
            )

    return evaluations_list


@pytest.fixture
def evaluation_scores(smiles_molw_config, evaluations):
    score_value = 3.0
    challenge = smiles_molw_config.challenge
    evaluation_score_type = challenge.scoretype_set.filter(
        level=models.ScoreType.Level.EVALUATION
    ).get()

    return [
        models.EvaluationScore.objects.create(
            submission_run=evaluation.submission_run,
            input_element=evaluation.input_element,
            score_type=evaluation_score_type,
            value=score_value,
        )
        for evaluation in evaluations
    ]
