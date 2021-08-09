import os.path

from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
from django.core.management.base import BaseCommand
from django.utils import timezone

from core import models

SAMPLE_INPUT_FILE = os.path.join(os.path.dirname(__file__), "data", "ChEBI_16716.mdl")

SAMPLE_OUTPUT_FILE = os.path.join(
    os.path.dirname(__file__), "data", "Conformer3D_CID_241.mdl"
)


def _create_challenge_inputs(challenge, file_based):
    if file_based:
        input_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=True,
            content_type=ContentType.objects.get_for_model(models.FileValue),
            key="molfile",
            description="MOL File",
        )

        output_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=False,
            content_type=ContentType.objects.get_for_model(models.FileValue),
            key="conformation",
            description="Conformation",
        )
    else:
        input_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=True,
            content_type=ContentType.objects.get_for_model(models.TextValue),
            key="smiles",
            description="SMILES",
        )

        output_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=False,
            content_type=ContentType.objects.get_for_model(models.FloatValue),
            key="molWeight",
            description="Molecular Weight",
        )

    for idx, name in enumerate(
        ("benzene", "phenylalanine", "octane", "toluene", "heptane", "hexane")
    ):
        elem = models.InputElement.objects.create(
            name=name, challenge=challenge, is_public=idx % 2
        )
        if file_based:
            value_object = models.FileValue.from_string(
                SAMPLE_INPUT_FILE, challenge=challenge
            )
            expected_value = models.FileValue.from_string(
                SAMPLE_OUTPUT_FILE, challenge=challenge
            )
        else:
            smiles = "c1ccccc1"
            value_object = models.TextValue.objects.create(
                challenge=challenge, value=smiles
            )
            expected_value = models.FloatValue.objects.create(
                challenge=challenge, value=72.0
            )

<<<<<<< HEAD
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
=======
>>>>>>> b9527f59ed08486be73de4b43aeb2be71e831c19
        models.InputValue.objects.create(
            input_element=elem, value_type=input_type, value_object=value_object
        )

        models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=output_type,
            value_object=expected_value,
        )


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("--name", default="Milo")
        parser.add_argument(
            "--files", action="store_true", help="File based sample data"
        )
        parser.add_argument(
            "--delete", action="store_true", help="Delete old challenges with this name"
        )

    def handle(self, *args, **options):
<<<<<<< HEAD
        print("handle() is being run")
=======
        if options["delete"]:
            self.stdout.write("Deleting old challenges", ending="")
            for challenge in models.Challenge.objects.filter(name=options["name"]):
                challenge.delete()
                self.stdout.write(".", ending="")
            self.stdout.write("Finished deleting")

>>>>>>> b9527f59ed08486be73de4b43aeb2be71e831c19
        User = get_user_model()

        (user, _) = User.objects.get_or_create(
            username="milo", email="braxton.robbason@gmail.com"
        )
        foo_url = "https://foo.com"
        challenge = models.Challenge.objects.create(
<<<<<<< HEAD
            name="Megs Test",
=======
            name=options["name"],
>>>>>>> b9527f59ed08486be73de4b43aeb2be71e831c19
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

        _create_challenge_inputs(challenge, options["files"])
