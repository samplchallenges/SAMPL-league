# pylint: skip-file
from re import S
from unittest.mock import Mock, patch

import pytest

from core import batching, models


@pytest.fixture
def batch_evaluation_factory(smiles_molw_config, input_elements):
    def batch_evaluation_maker(max_batch_size):
        challenge = smiles_molw_config.challenge
        challenge.max_batch_size = max_batch_size
        challenge.save()
        submission_run = smiles_molw_config.submission_run
        molw_type = smiles_molw_config.output_type

        ibg = challenge.current_batch_group()
        evaluations_list = []
        for input_batch in ibg.inputbatch_set.all():
            batch_evaluation = models.BatchEvaluation.objects.create(
                input_batch=input_batch, submission_run=submission_run
            )
            evaluations_list.append(batch_evaluation)
            value = "777.05"

            def fake_invert(_foo):
                for element in input_batch.elements():
                    yield element.id, value

            batcher = batching.BATCHERS["csv"]
            with patch.object(batcher, "invert", fake_invert):
                # "core.foo.FOO", fake_invert): #"core.batching.CSVBatch.invert", fake_invert):
                models.Prediction.load_batch_output(
                    challenge, batch_evaluation, molw_type, "foo"
                )
            batch_evaluation.batchup()
        return evaluations_list

    return batch_evaluation_maker


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


def _make_evaluation_scores(smiles_molw_config, input_elements):
    score_value = 3.0
    challenge = smiles_molw_config.challenge
    submission_run = smiles_molw_config.submission_run
    evaluation_score_type = challenge.scoretype_set.filter(
        level=models.ScoreType.Level.EVALUATION
    ).get()

    return [
        models.EvaluationScore.objects.create(
            submission_run=submission_run,
            input_element=input_element,
            score_type=evaluation_score_type,
            value=score_value,
        )
        for input_element in input_elements
    ]


@pytest.fixture
def evaluation_scores(smiles_molw_config, evaluations):
    return _make_evaluation_scores(
        smiles_molw_config, [evaluation.input_element for evaluation in evaluations]
    )


@pytest.fixture
def batch_evaluation_score_factory(smiles_molw_config):
    def batch_evaluation_score_maker(batch_evaluations):
        elements = [
            elem
            for batch_evaluation in batch_evaluations
            for elem in batch_evaluation.input_batch.elements()
        ]
        return _make_evaluation_scores(smiles_molw_config, elements)

    return batch_evaluation_score_maker
