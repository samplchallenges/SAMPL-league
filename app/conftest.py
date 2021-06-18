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
def challenge(db):
    empty_url = "http://github.com"
    start_at = datetime(2020, 1, 1, hour=1, tzinfo=timezone.utc)
    end_at = datetime(2020, 9, 1, hour=1, tzinfo=timezone.utc)
    challenge = models.Challenge(
        name="SAMPL1",
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
    return os.path.join("tests", "data")


@pytest.fixture
def benzene_from_mol(testing_data_path, challenge, molfile_type, molw_type, db):
    elem = models.InputElement.objects.create(
        name="benzene", challenge=challenge, is_public=True
    )
    molfile_name = os.path.join(testing_data_path, "ChEBI_16716.mdl")
    with open(molfile_name, "rb") as mol_fp:
        molfile_value = models.FileValue()
        molfile_value.value.save(molfile_name, File(mol_fp))
        molfile_value.save()

        # molfile_contents = mol_fp.read()
    # molfile_value = models.BlobValue.objects.create(value=molfile_contents)
    models.InputValue.objects.create(
        input_element=elem, value_type=molfile_type, value_object=molfile_value
    )
    float_value = models.FloatValue.objects.create(value=72.0)
    models.AnswerKey.objects.create(
        challenge=challenge,
        input_element=elem,
        value_type=molw_type,
        value_object=float_value,
    )
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
        smiles_value = models.TextValue.objects.create(value=smiles)
        models.InputValue.objects.create(
            input_element=elem, value_type=smiles_type, value_object=smiles_value
        )
        float_value = models.FloatValue.objects.create(value=72.0)
        models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=molw_type,
            value_object=float_value,
        )
        elems.append(elem)
    return elems
