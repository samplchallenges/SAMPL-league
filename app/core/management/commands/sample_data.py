import os

from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand
from django.utils import timezone

from core import models


def _create_challenge_inputs(challenge):
    elems = {}
    for idx, name in enumerate(
        ("benzene", "phenylalanine", "octane", "toluene", "heptane", "hexane")
    ):
        elems[name] = models.InputElement.objects.create(
            name=name, challenge=challenge, is_public=idx % 2
        )

    smiles_type = models.InputType.objects.create(
        challenge=challenge, key="SMILES", description="SMILES"
    )

    smiles = "c1cccc1"
    for elem in elems.values():
        models.InputValue.objects.create(
            input_element=elem, input_type=smiles_type, value=smiles
        )


class Command(BaseCommand):
    def handle(self, *args, **options):
        User = get_user_model()

        (user, _) = User.objects.get_or_create(
            username="milo", email="braxton.robbason@gmail.com"
        )
        foo_url = "https://foo.com"
        challenge = models.Challenge.objects.create(
            name="Milo",
            start_at=timezone.now(),
            end_at=timezone.now(),
            repo_url=foo_url,
            sample_data_url=foo_url,
            sample_score_reference_url=foo_url,
            secret_data_url=foo_url,
            secret_score_reference_url=foo_url,
            execution_options_json={"foo": "bar"},
        )

        container = models.Container.objects.create(
            name="crate",
            user=user,
            challenge=challenge,
            registry=foo_url,
            label="milo_crate",
        )

        scoring_container = models.Container.objects.create(
            name="scorer",
            user=user,
            challenge=challenge,
            registry=foo_url,
            label="milo_scorer",
        )

        scoremaker = models.ScoreMaker.objects.create(
            challenge=challenge, container=scoring_container
        )

        submission = models.Submission.objects.create(
            user=user,
            challenge=challenge,
            container=container,
            name="Milo's submission",
        )

        first_run = models.SubmissionRun.objects.create(
            submission=submission,
            digest="cafecafef00d",
            data_privacy_level=models.SubmissionRun._DataPrivacyLevel.PUBLIC,
            status=models.SubmissionRun._Status.SUCCESS,
        )

        _create_challenge_inputs(challenge)
