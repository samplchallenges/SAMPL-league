from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
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

    smiles_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=True,
        content_type=ContentType.objects.get_for_model(models.TextValue),
        key="smiles",
        description="SMILES",
    )

    molw_type = models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FloatValue),
        key="LogP",
        description="LogP Value",
    )

    smiles = "c1ccccc1"
    for elem in elems.values():
        smiles_value = models.TextValue.objects.create(
            challenge=challenge, value=smiles
        )
        models.InputValue.objects.create(
            input_element=elem, value_type=smiles_type, value_object=smiles_value
        )
        float_value = models.FloatValue.objects.create(challenge=challenge, value=72.0)
        models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=molw_type,
            value_object=float_value,
        )


class Command(BaseCommand):
    def handle(self, *args, **options):
        print("handle() is being run")
        User = get_user_model()

        (user, _) = User.objects.get_or_create(
            username="milo", email="braxton.robbason@gmail.com"
        )
        foo_url = "https://foo.com"
        challenge = models.Challenge.objects.create(
            name="Megs Test",
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
            registry="docker.io",
            label="osatom/logpcalc",
            tag="latest",
            # ghcr.io/robbason/calc-molwt:latest
            #registry="ghcr.io",
            #label="robbason/calc-coords",
            #tag="latest",
            # registry="docker.io",
            # label="mmh42/calc-molwt",
            # tag="0.1",
        )

        scoring_container = models.Container.objects.create(
            name="scorer",
            user=user,
            challenge=challenge,
            registry="docker.io",
            label="osatom/scoremore",
            tag="latest",
            #registry="ghcr.io",
            #label="robbason/score-coords",
            #tag="latest",
        )

        scoremaker = models.ScoreMaker.objects.create(
            challenge=challenge, container=scoring_container
        )

        models.ScoreType.objects.create(
            challenge=challenge, key="diff", level=models.ScoreType.Level.EVALUATION
        )

        models.ScoreType.objects.create(
            challenge=challenge, key="rmse", level=models.ScoreType.Level.SUBMISSION_RUN
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
            is_public=True,
            status=models.Status.SUCCESS,
        )

        _create_challenge_inputs(challenge)
