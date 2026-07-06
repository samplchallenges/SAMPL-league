import re
from functools import cached_property

import django.contrib.auth.models as auth_models
from django.conf import settings
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ValidationError
from django.db import models
from django.utils import timezone

from .. import configurator, filecache
from .infra_models import ContainerType, Timestamped


class NotFullyLoadedException(Exception):
    def __init__(self, object_description, errors):
        super().__init__()
        self.object_description = object_description
        self.errors = errors

    def __str__(self):
        errstr = "; ".join([f"{k}: {v}" for k, v in self.errors.items()])
        return f"Errors for {self.object_description}: {errstr}"


class Challenge(Timestamped):
    name = models.CharField(max_length=255, unique=True)
    start_at = models.DateTimeField()
    end_at = models.DateTimeField()
    repo_url = models.URLField()
    max_batch_size = models.PositiveIntegerField(
        help_text="0 to disable batching", default=0
    )

    __output_types_dict = None

    def __str__(self):
        return str(self.name)

    def is_active(self):
        return self.end_at > timezone.now() and self.start_at < timezone.now()

    def __load_output_types(self):
        output_types = self.valuetype_set.filter(is_input_flag=False)
        self.__output_types_dict = {
            output_type.key: output_type for output_type in output_types.all()
        }
        file_content_type = ContentType.objects.get_by_natural_key("core", "filevalue")
        self.__output_file_keys = {  # pylint:disable=attribute-defined-outside-init
            key
            for key, output_type in self.__output_types_dict.items()
            if output_type.content_type == file_content_type
        }

    def output_file_keys(self):
        if self.__output_types_dict is None:
            self.__load_output_types()
        return self.__output_file_keys

    def output_keys(self):
        if self.__output_types_dict is None:
            self.__load_output_types()
        return self.__output_types_dict.keys()

    def output_type(self, key):
        if self.__output_types_dict is None:
            self.__load_output_types()
        return self.__output_types_dict.get(key)

    def current_batch_group(self):
        if self.max_batch_size > 0:
            return self.inputbatchgroup_set.order_by("-created_at").first()
        return None

    @cached_property
    def score_types(self):
        score_types = {
            ScoreType.Level.EVALUATION: {},
            ScoreType.Level.SUBMISSION_RUN: {},
        }

        for score_type in self.scoretype_set.all():
            score_types[score_type.level][score_type.key] = score_type

        return score_types

    def fully_loaded(self):
        # TODO: slow in case of many input elements due to ORM queries being rerun
        errors = {}
        valid = True
        if not hasattr(self, "scoremaker"):
            valid = False
            errors["ScoreMaker"] = "Missing"

        if not self.score_types[ScoreType.Level.EVALUATION]:
            valid = False
            errors["Evaluation Score Types"] = "Missing"
        if not self.score_types[ScoreType.Level.SUBMISSION_RUN]:
            valid = False
            errors["Submission Run Score Types"] = "Missing"

        element_visibilities = set()
        for element in self.inputelement_set.all():
            element.fully_loaded()
            if not element.is_parent:
                element_visibilities.add(element.is_public)
        if element_visibilities != {False, True}:
            valid = False
            errors["Need both public and private elements"] = (
                f"Found {element_visibilities}"
            )
        if not valid:
            raise NotFullyLoadedException(str(self), errors)


class ScoreType(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    key = models.CharField(max_length=255)

    class Level(models.TextChoices):
        EVALUATION = "evaluation"
        SUBMISSION_RUN = "submission_run"

    level = models.CharField(
        max_length=255, choices=Level.choices, default=Level.EVALUATION
    )

    class Meta:
        unique_together = ["challenge", "key", "level"]

    def __str__(self):
        return f"{self.key}, {self.level}"


class Container(Timestamped):
    # Containers can be managed by admins or users, but keep with admin_managed
    #  since admins need to create scoring containers
    # TODO: perhaps put a base container class into infra_models and truly separate
    #  scoring containers from submission containers?
    name = models.CharField(max_length=255)
    user = models.ForeignKey(auth_models.User, on_delete=models.CASCADE)
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    container_type = models.CharField(
        max_length=255,
        choices=ContainerType.choices,
        help_text=configurator.CONTAINER_TYPE_DETAILS,
        null=True,
    )
    registry = models.CharField(max_length=255)
    label = models.CharField(max_length=255)
    tag = models.CharField(max_length=255, blank=True, null=True)
    digest = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self):
        return str(self.name)

    @property
    def uri(self):
        suffix = f":{self.tag}" if self.tag else ""
        # TODO: revert this change for local testing
        if settings.LOCAL_CONTAINERS:
            # return  "{self.label}{suffix}"
            return self.label
        return f"{self.registry}/{self.label}{suffix}"

    def custom_args(self):
        return {
            arg.key: arg.string_value for arg in self.args.all() if arg.string_value
        }

    def custom_file_args(self):
        return {
            arg.key: filecache.ensure_local_copy(arg.file_value)
            for arg in self.args.all()
            if arg.file_value
        }

    @property
    def local_save_path(self):
        if not settings.CONTAINER_FILES_ROOT:
            return

        def _sanitize_user_str(user_str):
            MAX_USER_STR_LEN = 20
            return re.sub(r"\W", "_", user_str[:MAX_USER_STR_LEN])

        clean_label = _sanitize_user_str(self.label)
        if self.tag:
            suffix = "_" + _sanitize_user_str(self.tag)
        else:
            suffix = ""
        return (
            settings.CONTAINER_FILES_ROOT
            / f"container_{self.id}"
            / f"{clean_label}{suffix}.sif"
        )


class ScoreMaker(Timestamped):
    """
    Each active challenge has exactly one scoring container.
    """

    challenge = models.OneToOneField(
        Challenge, on_delete=models.CASCADE, primary_key=True
    )
    container = models.ForeignKey(Container, on_delete=models.CASCADE)

    def __str__(self):
        return f"{self.challenge} {self.container}"


class InputElement(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    parent = models.ForeignKey(
        "self",
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        help_text="Inherit common values from parent element, eg. protein PDB file",
    )
    is_parent = models.BooleanField(
        default=False,
        help_text="Parent elements won't be evaluated. They are used to"
        " group common values",
    )
    name = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)

    class Meta:
        unique_together = ["challenge", "name"]

    def clean(self):
        if self.is_parent and self.parent is not None:
            raise ValidationError(
                "If the input element is a parent, it cannot also have a parent."
                " Only one level of parent-child relationship is supported"
            )

    def save_score(self, submission_run, score_type, value):
        self.evaluationscore_set.create(
            submission_run=submission_run, score_type=score_type, value=float(value)
        )
        print("EVAL SCORE COUNT:", self.evaluationscore_set.count())
        print("Evaluation scores:", self.evaluationscore_set.all())

    def loaded_input_keys(self):
        return {
            iv.value_type.key
            for iv in self.inputvalue_set.select_related("value_type").all()
        }

    def loaded_answer_keys(self):
        return {
            ak.value_type.key
            for ak in self.answerkey_set.select_related("value_type").all()
        }

    def fully_loaded(self):
        """
        Ensure all input values and answer keys are loaded for this input element
        """
        if self.is_parent:
            for element in self.inputelement_set.all():
                element.fully_loaded()
            return

        valid = True
        errors = dict()
        if self.parent:
            loaded_input_keys = self.parent.loaded_input_keys()
            loaded_answer_keys = self.parent.loaded_answer_keys()
        else:
            loaded_input_keys = set()
            loaded_answer_keys = set()
        loaded_input_keys.update(self.loaded_input_keys())
        expected_input_keys = {
            vt.key for vt in self.challenge.valuetype_set.filter(is_input_flag=True)
        }
        if loaded_input_keys != expected_input_keys:
            valid = False
            errors["Missing Input Keys"] = expected_input_keys - loaded_input_keys

        loaded_answer_keys.update(self.loaded_answer_keys())
        expected_answer_keys = self.challenge.output_keys()
        if loaded_answer_keys != expected_answer_keys:
            valid = False
            errors["Missing Answer Keys"] = expected_answer_keys - loaded_answer_keys

        if not valid:
            raise NotFullyLoadedException(str(self), errors)

    def __str__(self):
        return f"{self.name}, is public? {self.is_public}"


class ValueType(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    is_input_flag = models.BooleanField(choices=((True, "Input"), (False, "Output")))
    on_parent_flag = models.BooleanField(
        choices=((True, "On parent"), (False, "On child")),
        default=False,
        help_text="For input elements, should this value be set on the parent or on"
        " the child? e.g. protein PDB structure should only be set on the parent;"
        " ligand SMILES should only be set on the child",
    )
    key = models.CharField(max_length=255)
    batch_method = models.CharField(
        max_length=10,
        blank=True,
        default="csv",
        choices=(("", "None"), ("csv", "CSV"), ("sdf", "MOL or SDF")),
    )
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    description = models.TextField()

    class Meta:
        unique_together = ["challenge", "is_input_flag", "key"]

    def __str__(self):
        return self.key

    def type_input_note(self):
        if self.challenge.max_batch_size > 0 and not self.on_parent_flag:
            if self.batch_method == "csv":
                return "CSV file with columns ID, name, value"
            elif self.batch_method == "sdf":
                return "multi-molecule SDF with properties SAMPL_ID, SAMPL_NAME"
            else:
                return "Configuration error! Contact administrator."
        return self.content_type.name

    def type_output_note(self):
        if self.challenge.max_batch_size > 0:
            if self.batch_method == "csv":
                return "CSV file with columns ID, name, value, where the ID matches up with the values in the inputs"
            elif self.batch_method == "sdf":
                return "multi-molecule SDF with property SAMPL_ID that must match up with the inputs"
            else:
                return "Configuration error! Contact administrator."
        return self.content_type.name
