import os.path
import tempfile
from unittest.mock import Mock, patch

import pytest
from django.core.files import File

from core import filecache, models
from referee import scoring


def test_score_submission_run(
    smiles_molw_config,
    evaluations,
    evaluation_scores,
    container_engine,
):
    with patch("django.conf.settings.CONTAINER_ENGINE", container_engine):
        submission_run = smiles_molw_config.submission_run
        scoring.score_submission_run(submission_run)

        score_type = models.ScoreType.objects.get(
            challenge=smiles_molw_config.challenge, key="rmse"
        )

        assert models.SubmissionRunScore.objects.count() == 1
        submission_run_score = models.SubmissionRunScore.objects.get(
            score_type=score_type
        )
        assert submission_run_score.value == pytest.approx(3.0)
        assert models.EvaluationScore.objects.count() == 2

        evaluation = evaluations[0]
        assert evaluation.scores.count() == 1
        score = evaluation.scores.first()
        assert score.value == pytest.approx(3.0)


@pytest.mark.parametrize("max_batch_size", [1, 2])
def test_score_submission_run_batch(
    smiles_molw_config,
    batch_evaluation_factory,
    batch_evaluation_score_factory,
    container_engine,
    max_batch_size,
):
    batch_evaluations = batch_evaluation_factory(max_batch_size)
    scores = batch_evaluation_score_factory(batch_evaluations)
    with patch("django.conf.settings.CONTAINER_ENGINE", container_engine):
        submission_run = smiles_molw_config.submission_run
        scoring.score_submission_run(submission_run)

    score_type = models.ScoreType.objects.get(
        challenge=smiles_molw_config.challenge, key="rmse"
    )
    assert models.SubmissionRunScore.objects.count() == 1

    submission_run_score = models.SubmissionRunScore.objects.get(score_type=score_type)
    assert submission_run_score.value == pytest.approx(3.0)


def _save_file_arg(container, key, file_body):
    # custom_file_args has the file contents as value
    arg = models.ContainerArg(container=container, key=key)
    arg.save()
    with tempfile.TemporaryDirectory() as tmpdir:
        hellopath = os.path.join(tmpdir, "hello.txt")
        with open(hellopath, "w", encoding="utf-8") as fp:
            fp.write(file_body)
        with open(hellopath, "rb") as fp:
            arg.file_value.save(key, File(fp))
    arg.save()
    return arg


@pytest.mark.parametrize(
    ["custom_string_args", "custom_file_args"],
    [
        [{}, {}],
        [{"foo": "bar"}, {}],
        [{}, {"license": "Hello world"}],
    ],
)
def test_score_evaluation_args(
    smiles_molw_config,
    evaluations,
    custom_string_args,
    custom_file_args,
):
    submission_run = smiles_molw_config.submission_run
    challenge = submission_run.submission.challenge
    scoring_container = challenge.scoremaker.container
    evaluation_score_types = challenge.score_types[models.ScoreType.Level.EVALUATION]

    for k, v in custom_string_args.items():
        models.ContainerArg.objects.create(
            container=scoring_container, key=k, string_value=v
        )

    custom_file_args_full = {}
    for k, v in custom_file_args.items():
        file_arg = _save_file_arg(scoring_container, k, v)
        custom_file_args_full[k] = filecache.ensure_local_copy(file_arg.file_value)

    fake_run = Mock(return_value=[])

    evaluation = evaluations[0]
    assert evaluation.input_element.is_public

    with patch("referee.scoring.ever_given.wrapper") as mock_wrapper:
        mock_wrapper.run = fake_run

        scoring.score_element(
            scoring_container,
            evaluation,
            evaluation.submission_run,
            evaluation.input_element,
            evaluation_score_types,
        )

        expected_args = ("ghcr.io/robbason/score-coords:latest", "score-evaluation")
        expected_kwargs = {"molWeight_answerkey": 72.0, "molWeight_prediction": 97.08}
        expected_kwargs.update(custom_string_args)
        expected_file_kwargs = custom_file_args_full
        call_args = fake_run.call_args
        assert call_args is not None
        args = call_args.args
        assert args == expected_args
        kwargs = call_args.kwargs["kwargs"]
        assert kwargs == expected_kwargs
        file_kwargs = call_args.kwargs["file_kwargs"]
        assert file_kwargs == expected_file_kwargs


def test_score_submission_run_failure(
    smiles_molw_config,
    evaluations,
    evaluation_scores,
    container_engine,
):
    with patch("django.conf.settings.CONTAINER_ENGINE", container_engine):
        fake_run = Mock(return_value=[], side_effect=Exception("Planned Exception"))

        submission_run = smiles_molw_config.submission_run
        with patch("referee.scoring.ever_given.wrapper") as mock_wrapper:
            mock_wrapper.run = fake_run
            with pytest.raises(scoring.ScoringFailureException):
                scoring.score_submission_run(submission_run)


def test_score_evaluation_failure(smiles_molw_config, evaluations, container_engine):
    with patch("django.conf.settings.CONTAINER_ENGINE", container_engine):
        with patch("referee.scoring.ever_given.wrapper") as mock_wrapper:
            _evaluation_failure(smiles_molw_config, evaluations, mock_wrapper)


def _evaluation_failure(smiles_molw_config, evaluations, mock_wrapper):
    fake_run = Mock(return_value=[], side_effect=Exception("Planned Exception"))
    mock_wrapper.run = fake_run
    submission_run = smiles_molw_config.submission_run
    challenge = submission_run.submission.challenge
    scoring_container = challenge.scoremaker.container
    evaluation_score_types = challenge.score_types[models.ScoreType.Level.EVALUATION]
    evaluation = evaluations[0]

    with pytest.raises(scoring.ScoringFailureException):
        scoring.score_element(
            scoring_container,
            evaluation,
            evaluation.submission_run,
            evaluation.input_element,
            evaluation_score_types,
        )
