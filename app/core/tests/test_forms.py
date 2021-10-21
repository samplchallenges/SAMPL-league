# pylint: disable=unused-argument, unused-variable

import pytest
from datetime import timedelta

from core.forms import ContainerForm, SubmissionForm
from django.utils import timezone

from core import models


@pytest.mark.django_db
def test_create(challenge, user):
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

@pytest.mark.django_db
def test_expired_challenge(challenge, user):
    submission_form = SubmissionForm()
    assert not submission_form.is_valid()
    container_form = ContainerForm()
    assert not container_form.is_valid()

    challenge.start_at = timezone.now() - timedelta(hours=3)
    challenge.end_at = challenge.start_at + timedelta(hours=1)

    assert not challenge.is_active()
    assert challenge.end_at < timezone.now()


