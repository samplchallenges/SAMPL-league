from datetime import datetime
from datetime import timezone
import pytest

from django.contrib.auth import get_user_model

from core.models import Challenge, Container, Submission


@pytest.fixture
def user(db):
    User = get_user_model()
    user = User(username="hello")
    user.save()
    return user


@pytest.fixture
def challenge(db):
    empty_url = "http://github.com"
    start_at = datetime(2020, 1, 1, hour=1, tzinfo=timezone.utc)
    end_at = datetime(2020, 9, 1, hour=1, tzinfo=timezone.utc)
    challenge = Challenge(
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
    container = Container(
        name="Container1",
        user=user,
        challenge=challenge,
        registry="local",
        label="package1",
    )
    container.save()
    return container


@pytest.fixture
def draft_submission(container, db):
    submission = Submission(
        name="Draft Submission",
        user=container.user,
        container=container,
        challenge=container.challenge,
    )
    submission.save()
    return submission
