from datetime import datetime
from datetime import timezone
import pytest

from django.contrib.auth.models import User

from core.forms import SubmissionForm, ContainerForm
from core.models import Challenge

# run like so
# pytest --reuse-db tests/test_upload.py


@pytest.fixture
def user(db):
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


@pytest.mark.django_db
def test_create(challenge, user):
    submission_form = SubmissionForm()
    assert not submission_form.is_valid()
    container_form = ContainerForm()
    assert not container_form.is_valid()
    container_form = ContainerForm(
        data={
            "name": "My Container",
            "registry": "docker",
            "label": "foo",
            "challenge": challenge,
        }
    )
    print(container_form.errors)
    assert container_form.is_valid()
    container = container_form.save(commit=False)
    container.user = user
    container.save()
    submission_form = SubmissionForm(
        data={
            "name": "Test submission",
            "challenge": challenge,
            "container": container,
        }
    )
    print(submission_form.errors)
    assert submission_form.is_valid()

    submission = submission_form.save(commit=False)
    submission.user = user
    submission.save()
    assert submission.draft_mode
