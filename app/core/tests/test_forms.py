# pylint: disable=unused-argument, unused-variable

from unittest.mock import patch

import pytest

from core.forms import ContainerForm, SubmissionForm
from core.tests import mocktime


@patch("django.utils.timezone.now", mocktime.active)
@pytest.mark.django_db
def test_create(challenge, user):
    assert challenge.is_active()
    submission_form = SubmissionForm()
    assert not submission_form.is_valid()
    container_form = ContainerForm()
    assert not container_form.is_valid()

    container_form = ContainerForm(
        data={
            "container-name": "My Container",
            "container-challenge": challenge,
            "container-registry": "docker",
            "container-label": "foo",
        }
    )
    assert container_form.is_valid()
    container = container_form.save(commit=False)
    container.user = user
    container.save()
    submission_form = SubmissionForm(
        data={
            "submission-name": "Test submission",
        }
    )
    print(submission_form.errors)
    assert submission_form.is_valid()

    submission = submission_form.save(commit=False)
    submission.user = user
    submission.container = container
    submission.challenge = container.challenge
    submission.save()
    assert submission.draft_mode
