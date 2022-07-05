import datetime
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

SAMPLE_PDB_FILE = os.path.join(os.path.dirname(__file__), "data", "5qcr.pdb")


def _create_challenge_inputs(challenge, file_based, with_pdb):
    if with_pdb:
        pdb_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=True,
            on_parent_flag=True,
            content_type=ContentType.objects.get_for_model(models.FileValue),
            key="protein_pdb",
            description="Protein structure (shared among input elements)",
        )
        parent_elem = models.InputElement.objects.create(
            name="Protein", challenge=challenge, is_public=True, is_parent=True
        )
        value_object = models.FileValue.from_string(
            SAMPLE_PDB_FILE, challenge=challenge, input_element=parent_elem
        )
        value_object.save()
        models.InputValue.objects.create(
            input_element=parent_elem,
            value_type=pdb_type,
            value_object=value_object,
        )
    else:
        parent_elem = None

    if file_based:
        input_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=True,
            content_type=ContentType.objects.get_for_model(models.FileValue),
            key="molfile",
            description="SD File",
            batch_method="sdf",
        )

        output_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=False,
            content_type=ContentType.objects.get_for_model(models.FileValue),
            key="conformation",
            description="Conformation",
            batch_method="sdf",
        )
    else:
        input_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=True,
            content_type=ContentType.objects.get_for_model(models.TextValue),
            key="smiles",
            description="SMILES",
            batch_method="csv",
        )

        output_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=False,
            content_type=ContentType.objects.get_for_model(models.FloatValue),
            key="molWeight",
            description="Molecular Weight",
            batch_method="csv",
        )

    for idx, name in enumerate(
        ("benzene", "phenylalanine", "octane", "toluene", "heptane", "hexane")
    ):
        elem = models.InputElement.objects.create(
            name=name, challenge=challenge, parent=parent_elem, is_public=idx % 2
        )
        if file_based:
            value_object = models.FileValue.from_string(
                SAMPLE_INPUT_FILE, challenge=challenge, input_element=elem
            )
            value_object.save()
            expected_value = models.FileValue.from_string(
                SAMPLE_OUTPUT_FILE, challenge=challenge, input_element=elem
            )
            expected_value.save()
        else:
            smiles = "c1ccccc1"
            value_object = models.TextValue.objects.create(
                challenge=challenge, value=smiles, input_element=elem
            )
            expected_value = models.FloatValue.objects.create(
                challenge=challenge, value=72.0, input_element=elem
            )

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
        parser.add_argument(
            "--withpdb", action="store_true", help="Show how to pass PDB file"
        )

    def handle(self, *args, **options):  # pylint:disable=unused-argument
        if options["delete"]:
            self.stdout.write("Deleting old challenges", ending="")
            for challenge in models.Challenge.objects.filter(name=options["name"]):
                challenge.delete()
                self.stdout.write(".", ending="")
            self.stdout.write("Finished deleting")

        User = get_user_model()

        (user, _) = User.objects.get_or_create(
            username="milo", email="braxton.robbason@gmail.com"
        )
        foo_url = "https://foo.com"
        challenge, _ = models.Challenge.objects.get_or_create(
            name=options["name"],
            defaults={
                "start_at": timezone.now(),
                "end_at": timezone.now() + datetime.timedelta(days=100),
                "repo_url": foo_url,
            },
        )

        container = models.Container.objects.create(
            name="crate",
            user=user,
            challenge=challenge,
            container_type="docker",
            registry="ghcr.io",
            label="robbason/calc-coords",
            tag="latest",
        )

        scoring_container = models.Container.objects.create(
            name="scorer",
            user=user,
            challenge=challenge,
            container_type="docker",
            registry="ghcr.io",
            label="megosato/score-coords",
            tag="latest",
        )

        models.ScoreMaker.objects.create(
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

        models.SubmissionRun.objects.create(
            submission=submission,
            digest="cafecafef00d",
            is_public=True,
            status=models.Status.SUCCESS,
        )

        _create_challenge_inputs(challenge, options["files"], options["withpdb"])
