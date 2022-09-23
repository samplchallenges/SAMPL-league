import csv
from pathlib import Path

from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
from django.core.management.base import BaseCommand
from django.db import transaction
from django.utils.timezone import make_aware
from ruamel.yaml import YAML

from core import models

STORAGE_TO_CLS = {
    "file": models.FileValue,
    "string": models.TextValue,
    "number": models.FloatValue,
}


def _content_type(storage):
    cls = STORAGE_TO_CLS[storage]
    return ContentType.objects.get_for_model(cls)


def _create_input_types(challenge, itconfigs):
    return _create_types(challenge, itconfigs, is_input_flag=True)


def _create_output_types(challenge, itconfigs):
    return _create_types(challenge, itconfigs, is_input_flag=False)


def _create_types(challenge, itconfigs, **kwargs):
    types = {}
    for itconfig in itconfigs:
        content_type = _content_type(itconfig["storage"])
        del itconfig["storage"]
        types[itconfig["key"]] = models.ValueType.objects.create(
            challenge=challenge, content_type=content_type, **itconfig, **kwargs
        )
    return types


def _create_values(
    challenge, config_path, elem, input_types, output_types, values_dict
):
    for key, value in values_dict.items():
        input_type = input_types.get(key)
        if input_type:
            _create_input_value(challenge, config_path, input_type, elem, value)
        else:
            output_type = output_types.get(key)
            _create_output_value(challenge, config_path, output_type, elem, value)


def _create_value(challenge, config_path, value_type, input_element, value):
    storage_cls = value_type.content_type.model_class()
    # TODO: to support multi-mol SDF inputs, need to:
    # look up input element name in attributes of file when loading
    # OR allow SDF input instead of CSV in this loop - prob more elegant!
    if storage_cls == models.FileValue:
        # Allow relative paths in config yaml
        # raise Exception(f"{config_path}, {config_path.parent}, {value}, {config_path.parent / value}")
        value = config_path.parent / value
    value_object = storage_cls.from_string(
        value, challenge=challenge, input_element=input_element
    )
    value_object.save()
    return value_object


def _create_input_value(challenge, config_path, value_type, input_element, value):
    value_object = _create_value(
        challenge, config_path, value_type, input_element, value
    )
    models.InputValue.objects.create(
        input_element=input_element,
        value_type=value_type,
        value_object=value_object,
    )


def _create_output_value(challenge, config_path, value_type, input_element, value):
    value_object = _create_value(
        challenge, config_path, value_type, input_element, value
    )
    models.AnswerKey.objects.create(
        challenge=challenge,
        input_element=input_element,
        value_type=value_type,
        value_object=value_object,
    )


def _create_elements_from_parents(
    challenge, config_path, elconfigs, input_types, output_types
):
    for elconfig in elconfigs:
        parent_element = models.InputElement.objects.create(
            challenge=challenge,
            is_parent=True,
            is_public=elconfig["public"],
            name=elconfig["name"],
        )
        for key, value in elconfig["inputs"].items():
            value_type = input_types[key]
            _create_input_value(
                challenge, config_path, value_type, parent_element, value
            )
        for key, value in elconfig["outputs"].items():
            value_type = output_types[key]
            _create_output_value(
                challenge, config_path, value_type, parent_element, value
            )
        _create_elements_from_files(
            challenge,
            config_path,
            elconfig["element_files"],
            input_types,
            output_types,
            parent_element=parent_element,
        )


def _create_elements_from_files(
    challenge, config_path, elconfigs, input_types, output_types, parent_element=None
):
    for element_file in elconfigs:
        is_public = element_file["public"]
        _create_elements_from_file(
            challenge,
            config_path,
            element_file["path"],
            input_types,
            output_types,
            parent_element=parent_element,
            is_public=is_public,
        )


def _create_elements_from_file(
    challenge,
    config_path,
    filename,
    input_types,
    output_types,
    parent_element=None,
    is_public=False,
):
    is_validated = False
    with open(config_path.parent / filename, encoding="utf8") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            if not is_validated:
                in_keys = [
                    input_type.key
                    for input_type in input_types.values()
                    if not input_type.on_parent_flag
                ]
                out_keys = [
                    output_type.key
                    for output_type in output_types.values()
                    if not output_type.on_parent_flag
                ]
                expected_keys = {"name", *in_keys, *out_keys}
                if row.keys() != expected_keys:
                    raise ValueError(
                        f"CSV issue with {filename}; expected keys {expected_keys} but found {row.keys()}"
                    )

            elem = models.InputElement.objects.create(
                name=row["name"],
                challenge=challenge,
                parent=parent_element,
                is_public=is_public,
            )
            del row["name"]
            _create_values(challenge, config_path, elem, input_types, output_types, row)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("filename")
        parser.add_argument(
            "--delete", action="store_true", help="Delete old challenges with this name"
        )
        parser.add_argument("--owner", default="admin", help="Scoring container owner")

    @transaction.atomic
    def handle(self, *args, **options):  # pylint:disable=unused-argument
        config_path = Path(options["filename"])
        yaml = YAML(typ="safe")  # default, if not specfied, is 'rt' (round-trip)
        with open(config_path, encoding="utf8") as fp:
            config_dict = yaml.load(fp)
        challenge_name = config_dict["name"]
        if options["delete"]:
            self.stdout.write("Deleting old challenges", ending="")
            for challenge in models.Challenge.objects.filter(name=challenge_name):
                challenge.delete()
                self.stdout.write(".", ending="")
            self.stdout.write(" Finished deleting")
        if models.Challenge.objects.filter(name=challenge_name).exists():
            raise Exception(
                f"Challenge named {challenge_name} already exists, run with --delete"
            )
        User = get_user_model()

        (user, _) = User.objects.get_or_create(
            username=options["owner"], email="braxton.robbason@gmail.com"
        )
        challenge_defaults = {
            k: v for k, v in config_dict.items() if k in ("repo_url", "max_batch_size")
        }
        challenge_dates = {
            k: make_aware(v)
            for k, v in config_dict.items()
            if k in ("start_at", "end_at")
        }
        challenge_defaults.update(challenge_dates)
        challenge, _ = models.Challenge.objects.get_or_create(
            name=challenge_name, defaults=challenge_defaults
        )

        scoring_container = models.Container.objects.create(
            user=user,
            challenge=challenge,
            name=f"Scoring Container for {challenge.name}",
            **config_dict["scoring"]["container"],
        )

        models.ScoreMaker.objects.create(
            challenge=challenge, container=scoring_container
        )

        for score_key in config_dict["scoring"]["element_types"]:
            models.ScoreType.objects.create(
                challenge=challenge,
                key=score_key,
                level=models.ScoreType.Level.EVALUATION,
            )
        for score_key in config_dict["scoring"]["run_types"]:
            models.ScoreType.objects.create(
                challenge=challenge,
                key=score_key,
                level=models.ScoreType.Level.SUBMISSION_RUN,
            )

        input_types = _create_input_types(challenge, config_dict["input_types"])
        output_types = _create_output_types(challenge, config_dict["output_types"])
        if "parent_elements" in config_dict:
            _create_elements_from_parents(
                challenge,
                config_path,
                config_dict["parent_elements"],
                input_types,
                output_types,
            )
        else:
            _create_elements_from_files(
                challenge,
                config_path,
                config_dict["element_files"],
                input_types,
                output_types,
            )
