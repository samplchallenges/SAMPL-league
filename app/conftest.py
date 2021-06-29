import os
from datetime import datetime, timezone

import dask.distributed as dd
import pytest
from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
from django.core.files import File

from core import models


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
            sample_data_url=empty_url,
            sample_score_reference_url=empty_url,
            secret_data_url=empty_url,
            secret_score_reference_url=empty_url,
            execution_options_json={},
        )
        challenge.save()
        return challenge
    return maker


@pytest.fixture
def container(challenge, user, db):
    return models.Container.objects.create(
        name="Container1",
        user=user,
        challenge=challenge,
        registry="ghcr.io",
        label="robbason/calc-molwt",
        tag="latest",
    )


@pytest.fixture
def scoring_container(challenge, user, db):
    return models.Container.objects.create(
        name="subtraction container",
        user=user,
        challenge=challenge,
        registry="docker.io",
        label="mmh42/calc-subtract",
        tag="0.1",
    )


@pytest.fixture
def score_maker(challenge, scoring_container, db):
    return models.ScoreMaker.objects.create(
        challenge=challenge, container=scoring_container
    )


@pytest.fixture
def score_types(challenge, db):
    return [
        models.ScoreType.objects.create(challenge=challenge, key=key, level=level)
        for key, level in (
            ("diff", models.ScoreType.Level.EVALUATION),
            ("rmse", models.ScoreType.Level.SUBMISSION_RUN),
        )
    ]


@pytest.fixture
def draft_submission(container, db):
    return models.Submission.objects.create(
        name="Draft Submission",
        user=container.user,
        container=container,
        challenge=container.challenge,
    )


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
        with open(file_path, "rb") as fp:
            file_value = models.FileValue()
            file_value.value.save(file_name, File(fp))
            file_value.save()
        models.InputValue.objects.create(
            input_element=elem, value_type=file_type, value_object=file_value
        )
        return elem
    return elem_maker


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
        with open(file_path, "rb") as fp:
            file_value = models.FileValue()
            file_value.value.save(file_name, File(fp))
            file_value.save()
        answer_key = models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=value_type,
            value_object=file_value,
        )
        return answer_key
    return fak_maker


@pytest.fixture
def benzene_from_mol(challenge, molfile_type, molw_type, elem_factory,
                     float_answer_key_factory, db):
    elem = elem_factory(challenge, molfile_type, "benzene", "ChEBI_16716.mdl")
    answer_key = float_answer_key_factory(challenge, elem, molw_type, 72.0)
    return elem


@pytest.fixture
def smiles_type(challenge, db):
    return models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=True,
        content_type=ContentType.objects.get_for_model(models.TextValue),
        key="SMILES",
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
def input_elements(challenge, smiles_type, molw_type, db):
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
