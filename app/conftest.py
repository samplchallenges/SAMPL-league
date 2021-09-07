import os
from collections import namedtuple
from datetime import datetime, timezone

import dask.distributed as dd
import pytest
from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
from django.core.files import File

from core import models

# Collection tuple to simplify test setup
TestConfig = namedtuple(
    "TestConfig", ["challenge", "input_type", "output_type", "submission_run"]
)


@pytest.fixture(scope="session")
def dask_client():
    # We need at least 4 workers
    with dd.LocalCluster(n_workers=4, preload=("daskworkerinit_tst.py",)) as cluster:
        yield dd.Client(cluster)


@pytest.fixture
def user(db):
    User = get_user_model()
    user, _ = User.objects.get_or_create(username="hello")
    return user


@pytest.fixture
def challenge(challenge_factory, db):
    return challenge_factory("SAMPL1")


@pytest.fixture
def challenge_factory(db):
    def maker(name):
        empty_url = "http://github.com"
        start_at = datetime(2020, 1, 1, hour=1, tzinfo=timezone.utc)
        end_at = datetime(2020, 9, 1, hour=1, tzinfo=timezone.utc)
        challenge = models.Challenge(
            name=name,
            start_at=start_at,
            end_at=end_at,
            repo_url=empty_url,
        )
        challenge.save()
        return challenge

    return maker


@pytest.fixture
def config_factory(challenge_factory, container_factory, db):
    """
    Create a challenge and related objects for testing
    """

    def maker(label, score_label, input_key, input_model, output_key, output_model):
        challenge = challenge_factory("standard")
        scoring_container = container_factory(challenge, score_label, tag="latest")
        models.ScoreMaker.objects.create(
            challenge=challenge, container=scoring_container
        )
        for key, level in (
            ("diff", models.ScoreType.Level.EVALUATION),
            ("rmse", models.ScoreType.Level.SUBMISSION_RUN),
        ):
            models.ScoreType.objects.create(challenge=challenge, key=key, level=level)

        input_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=True,
            content_type=ContentType.objects.get_for_model(input_model),
            key=input_key,
            description="N/A",
        )
        output_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=False,
            content_type=ContentType.objects.get_for_model(output_model),
            key=output_key,
            description="N/A",
        )
        container = container_factory(challenge, label, tag="latest")
        submission = models.Submission.objects.create(
            name="Draft Submission",
            user=container.user,
            container=container,
            challenge=challenge,
        )
        submission_run = models.SubmissionRun.objects.create(
            submission=submission,
            digest="cafef00d",
            is_public=True,
            status=models.Status.PENDING,
        )
        return TestConfig(challenge, input_type, output_type, submission_run)

    return maker


@pytest.fixture
def container_factory(user, db):
    def container_maker(challenge, label, tag):
        return models.Container.objects.create(
            name="Container1",
            user=user,
            challenge=challenge,
            registry="ghcr.io",
            label=label,
            tag=tag,
        )

    return container_maker


@pytest.fixture
def container(container_factory, challenge):
    return container_factory(
        challenge,
        label="robbason/calc-molwt",
        tag="latest",
    )


@pytest.fixture
def scoring_container(challenge, user, db):
    return models.Container.objects.create(
        name="subtraction container",
        user=user,
        challenge=challenge,
        registry="ghcr.io",
        label="robbason/score-coords",
        tag="latest",
    )


@pytest.fixture
def score_maker(challenge, scoring_container, db):
    return models.ScoreMaker.objects.create(
        challenge=challenge, container=scoring_container
    )


@pytest.fixture
def draft_submission(submission_factory, container, db):
    return submission_factory(container)


@pytest.fixture
def molfile_type(challenge, db):
    return models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=True,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="molfile",
        description="MOL file",
    )


@pytest.fixture
def testing_data_path():
    return os.path.join(os.path.dirname(__file__), "tests", "data")


@pytest.fixture
def elem_factory(testing_data_path, db):
    def elem_maker(challenge, file_type, name, file_name):
        elem = models.InputElement.objects.create(
            name=name, challenge=challenge, is_public=True
        )
        file_path = os.path.join(testing_data_path, file_name)
        file_value = models.FileValue.from_string(file_path, challenge=challenge)
        file_value.save()
        models.InputValue.objects.create(
            input_element=elem, value_type=file_type, value_object=file_value
        )
        return elem

    return elem_maker


@pytest.fixture
def submission_factory(db):
    def submission_maker(container):
        return models.Submission.objects.create(
            name="Draft Submission",
            user=container.user,
            container=container,
            challenge=container.challenge,
        )

    return submission_maker


@pytest.fixture
def submission_run_factory(db):
    def srun_maker(submission):
        return models.SubmissionRun.objects.create(
            submission=submission,
            digest="cafef00d",
            is_public=True,
            status=models.Status.PENDING,
        )

    return srun_maker


@pytest.fixture
def float_answer_key_factory(db):
    def fak_maker(challenge, elem, value_type, value):
        float_value = models.FloatValue.from_string("72.0", challenge=challenge)
        float_value.save()
        answer_key = models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=value_type,
            value_object=float_value,
        )
        return answer_key

    return fak_maker


@pytest.fixture
def file_answer_key_factory(testing_data_path, db):
    def fak_maker(challenge, elem, value_type, file_name):
        file_path = os.path.join(testing_data_path, file_name)
        file_value = models.FileValue.from_string(file_path, challenge=challenge)
        answer_key = models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=value_type,
            value_object=file_value,
        )
        return answer_key

    return fak_maker


@pytest.fixture
def molfile_molw_config(config_factory):
    return config_factory(
        "robbason/calc-molwt",
        "robbason/score-coords",
        "molfile",
        models.FileValue,
        "molWeight",
        models.FloatValue,
    )


@pytest.fixture
def smiles_molw_config(config_factory):
    return config_factory(
        "robbason/calc-molwt",
        "robbason/score-coords",
        "smiles",
        models.TextValue,
        "molWeight",
        models.FloatValue,
    )


@pytest.fixture
def benzene_from_mol(molfile_molw_config, elem_factory, float_answer_key_factory, db):
    elem = elem_factory(
        molfile_molw_config.challenge,
        molfile_molw_config.input_type,
        "benzene",
        "ChEBI_16716.mdl",
    )
    answer_key = float_answer_key_factory(
        molfile_molw_config.challenge, elem, molfile_molw_config.output_type, 72.0
    )
    return elem


@pytest.fixture
def smiles_type(challenge, db):
    return models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=True,
        content_type=ContentType.objects.get_for_model(models.TextValue),
        key="smiles",
        description="SMILES",
    )


@pytest.fixture
def molw_type(challenge, db):
    return models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FloatValue),
        key="molWeight",
        description="Molecular Weight",
    )


@pytest.fixture
def input_elements(smiles_molw_config, db):
    challenge = smiles_molw_config.challenge
    smiles_type = smiles_molw_config.input_type
    molw_type = smiles_molw_config.output_type
    elems = []
    for idx, (name, smiles) in enumerate(
        [
            ("benzene", "c1ccccc1"),
            ("octane", "C" * 8),
            ("hexane", "C" * 6),
            ("phenol", "c1ccccc1O"),
        ]
    ):
        elem = models.InputElement.objects.create(
            name=name, challenge=challenge, is_public=idx % 2
        )
        smiles_value = models.TextValue.from_string(smiles, challenge=challenge)
        smiles_value.save()
        models.InputValue.objects.create(
            input_element=elem, value_type=smiles_type, value_object=smiles_value
        )
        float_value = models.FloatValue.from_string("72.0", challenge=challenge)
        float_value.save()
        models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=molw_type,
            value_object=float_value,
        )
        elems.append(elem)
    return elems
