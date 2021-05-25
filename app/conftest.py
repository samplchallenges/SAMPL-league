from datetime import datetime, timezone

import dask.distributed as dd
import pytest
from django.contrib.auth import get_user_model

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
        registry="docker.io",
        label="mmh42/calc-molwt",
        tag="0.1",
    )


@pytest.fixture
def draft_submission(container, db):
    return models.Submission.objects.create(
        name="Draft Submission",
        user=container.user,
        container=container,
        challenge=container.challenge,
    )


@pytest.fixture
def smiles_input_type(challenge, db):
    return models.InputType.objects.create(
        challenge=challenge, key="SMILES", description="SMILES"
    )


@pytest.fixture
def input_elements(challenge, smiles_input_type, db):
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
        models.InputValue.objects.create(
            input_element=elem, input_type=smiles_input_type, value=smiles
        )
        float_value = models.FloatValue.objects.create(value=72.0)
        models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            key="molWeight",
            value_object=float_value,
        )
        elems.append(elem)
    return elems
