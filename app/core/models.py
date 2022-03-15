import logging
import os.path
import time
from collections import namedtuple
from functools import cached_property

import django.contrib.auth.models as auth_models
from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ValidationError
from django.core.files import File
from django.db import models, transaction
from django.db.models.functions import Concat
from django.urls import reverse
from django.utils import timezone
from django.utils.translation import gettext_lazy as _
from ever_given.utils import LogHandlerBase

from . import configurator, filecache

logger = logging.getLogger(__name__)


Completion = namedtuple("Completion", ["completed", "not_completed", "completed_frac"])


class Status(models.TextChoices):
    FAILURE = "FAILURE"
    SUCCESS = "SUCCESS"
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    CANCELLED = "CANCELLED"
    CANCEL_PENDING = "CANCEL_PENDING"


class StatusMixin:
    def is_finished(self):
        return self.status in (Status.SUCCESS, Status.FAILURE, Status.CANCELLED)


class Timestamped(models.Model):
    created_at = models.DateTimeField(default=timezone.now)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        abstract = True
        get_latest_by = ("created_at",)


def _timestamped_log(log):
    return " ".join([time.strftime("[%c %Z]", time.gmtime()), log])


class Logged(Timestamped, StatusMixin):
    log_stdout = models.TextField(blank=True)
    log_stderr = models.TextField(blank=True)

    log_handler = None

    class Meta:
        abstract = True

    class LogHandler(LogHandlerBase):
        def __init__(self, instance):
            self.instance_id = instance.pk
            self.cls = type(instance)

        @property
        def _update(self):
            return self.cls.objects.filter(pk=self.instance_id).update

        def handle_stdout(self, log):
            self._update(
                log_stdout=Concat(
                    models.F("log_stdout"), models.Value(_timestamped_log(log))
                )
            )

        def handle_stderr(self, log):
            self._update(
                log_stderr=Concat(
                    models.F("log_stderr"), models.Value(_timestamped_log(log))
                )
            )

    def append(self, stdout=None, stderr=None):
        update_fields = []
        if stdout is not None:
            self.log_stdout = Concat(models.F("log_stdout"), models.Value(stdout))
            update_fields.append("log_stdout")
        if stderr is not None:
            self.log_stderr = Concat(models.F("log_stderr"), models.Value(stderr))
            update_fields.append("log_stderr")

        if update_fields:
            self.save(update_fields=update_fields)


class Challenge(Timestamped):
    name = models.CharField(max_length=255, unique=True)
    start_at = models.DateTimeField()
    end_at = models.DateTimeField()
    repo_url = models.URLField()

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
        file_content_type = ContentType.objects.get_for_model(FileValue)
        self.__output_file_keys = {  # pylint:disable=attribute-defined-outside-init
            key
            for key, output_type in self.__output_types_dict.items()
            if output_type.content_type == file_content_type
        }

    def output_file_keys(self):
        if self.__output_types_dict is None:
            self.__load_output_types()
        return self.__output_file_keys

    def output_type(self, key):
        if self.__output_types_dict is None:
            self.__load_output_types()
        return self.__output_types_dict.get(key)

    @cached_property
    def score_types(self):
        score_types = {
            ScoreType.Level.EVALUATION: {},
            ScoreType.Level.SUBMISSION_RUN: {},
        }

        for score_type in self.scoretype_set.all():
            score_types[score_type.level][score_type.key] = score_type

        return score_types


class Container(Timestamped):
    name = models.CharField(max_length=255)
    user = models.ForeignKey(auth_models.User, on_delete=models.CASCADE)
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    registry = models.CharField(max_length=255)
    label = models.CharField(max_length=255)
    tag = models.CharField(max_length=255, blank=True, null=True)
    digest = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self):
        return str(self.name)

    @property
    def uri(self):
        suffix = f":{self.tag}" if self.tag else ""
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


class Submission(Timestamped):
    user = models.ForeignKey(auth_models.User, on_delete=models.CASCADE)
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    container = models.ForeignKey(Container, on_delete=models.CASCADE)
    draft_mode = models.BooleanField(
        default=False, help_text=configurator.DRAFT_MODE_DETAILS
    )
    # From user:
    name = models.CharField(max_length=40, help_text=configurator.NAME_DETAILS)
    url = models.URLField(blank=True, null=True)
    compute_time = models.TextField(
        help_text=configurator.COMPUTE_TIME_DETAILS, blank=True, null=True
    )
    computing_and_hardware = models.TextField(
        help_text=configurator.COMPUTING_AND_HARDWARE_DETAILS, blank=True, null=True
    )
    software = models.TextField(
        help_text=configurator.SOFTWARE_DETAILS, blank=True, null=True
    )
    category = models.CharField(
        max_length=255,
        choices=configurator.CATEGORY_CHOICES,
        help_text=configurator.CATEGORY_DETAILS,
        blank=True,
        null=True,
    )
    method = models.TextField(
        help_text=configurator.METHOD_DETAILS, blank=True, null=True
    )
    ranked = models.BooleanField(default=True, help_text=configurator.RANKED_DETAILS)

    notes = models.TextField(
        help_text=configurator.NOTES_DETAILS, blank=True, null=True
    )

    NONDRAFT_FIELD_NAMES = [
        "ranked",
        "category",
        "url",
        "compute_time",
        "computing_and_hardware",
        "software",
        "method",
    ]
    DETAIL_FIELD_NAMES = [*NONDRAFT_FIELD_NAMES, "notes"]

    def __str__(self):
        return f"{self.user}: {self.challenge}: {self.name}"

    def get_absolute_url(self):
        return reverse("submission", kwargs={"pk": self.pk})

    @property
    def missing_fields(self):
        return [
            field_name
            for field_name in self.NONDRAFT_FIELD_NAMES
            if getattr(self, field_name) in (None, "")
        ]

    @property
    def details(self):
        return [
            (field_name, getattr(self, field_name))
            for field_name in self.DETAIL_FIELD_NAMES
            if getattr(self, field_name) not in (None, "")
        ]

    def clean(self):
        super().clean()
        if not self.draft_mode:
            # All fields
            if self.missing_fields:
                self.draft_mode = True

    def clone(self):
        """
        Return a new submission using this as a template.
        Set draft mode and clear 'Name' field.
        """
        self.id = None  # pylint:disable=attribute-defined-outside-init
        self.pk = None  # pylint:disable=attribute-defined-outside-init
        self._state.adding = True
        return self

    def last_run(self, is_public):
        # TODO: this is inefficient the way it's used in submission_list
        return (
            self.submissionrun_set.filter(is_public=is_public)
            .order_by("-updated_at")
            .first()
        )

    def last_public_run(self):
        return self.last_run(is_public=True)

    def last_private_run(self):
        return self.last_run(is_public=False)

    def create_run(self, *, is_public):
        digest = self.container.digest
        if digest is None:
            digest = "nodigest"
        submission_run = self.submissionrun_set.create(
            digest=digest,
            is_public=is_public,
            status=Status.PENDING,
        )
        for element_id in self.challenge.inputelement_set.filter(
            is_public=is_public
        ).values_list("id", flat=True):
            submission_run.evaluation_set.create(input_element_id=element_id)

        return submission_run


def _container_file_location(instance, filename):
    return os.path.join(
        "container_args",
        str(instance.container.user_id),
        str(instance.container.id),
        instance.key,
        filename,
    )


class ContainerArg(Timestamped):
    container = models.ForeignKey(
        Container, on_delete=models.CASCADE, related_name="args"
    )
    key = models.SlugField(db_index=False)
    string_value = models.TextField(blank=True, null=True)
    file_value = models.FileField(
        blank=True, null=True, upload_to=_container_file_location
    )

    class Meta:
        unique_together = ["container", "key"]

    def clean(self):
        if not (self.string_value or self.file_value) or (
            self.string_value and self.file_value
        ):
            error = "Exactly one of string_value or file_value must be set"
            raise ValidationError(error)

    def filename(self):
        return os.path.basename(self.file_value.name)


class SubmissionRun(Logged):
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    digest = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)
    status = models.CharField(max_length=25, choices=Status.choices)

    def __str__(self):
        return f"{self.submission}:{self.digest}, status {self.status}"

    def check_cancel_requested(self):
        status = SubmissionRun.objects.values_list("status", flat=True).get(pk=self.id)
        return status == Status.CANCEL_PENDING

    def mark_for_cancel(self):
        with transaction.atomic():
            SubmissionRun.objects.filter(pk=self.id).update(
                status=Status.CANCEL_PENDING
            )
            Evaluation.objects.filter(
                submission_run_id=self.id, status=Status.PENDING
            ).update(status=Status.CANCELLED)
            Evaluation.objects.filter(
                submission_run_id=self.id, status=Status.RUNNING
            ).update(status=Status.CANCEL_PENDING)

    def completion(self):
        completed = self.evaluation_set.filter(
            status__in=(Status.SUCCESS, Status.FAILURE)
        ).count()
        num_element_ids = self.submission.challenge.inputelement_set.filter(
            is_public=self.is_public
        ).count()
        completed_frac = completed / num_element_ids if num_element_ids else 0
        return Completion(completed, num_element_ids, completed_frac)


class InputElement(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    name = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)

    class Meta:
        unique_together = ["challenge", "name"]

    @classmethod
    def all_element_values(cls, input_values):
        file_values = {}
        values = {}

        for input_value in input_values:
            value_type = input_value.value_type
            key = value_type.key
            content_type = value_type.content_type
            if content_type.model_class() == FileValue:
                file_values[key] = filecache.ensure_local_copy(input_value.value)
            else:
                values[key] = input_value.value
        return values, file_values

    def all_values(self):
        """
        Returns a pair of key: value dicts, where the first dict is regular values
        and the second is file values
        """
        input_values = self.inputvalue_set.select_related(
            "value_type", "value_type__content_type"
        ).all()
        return self.all_element_values(input_values)

    def __str__(self):
        return f"{self.name}, is public? {self.is_public}"


class ValueType(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    is_input_flag = models.BooleanField(choices=((True, "Input"), (False, "Output")))
    key = models.CharField(max_length=255)
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    description = models.TextField()

    class Meta:
        unique_together = ["challenge", "is_input_flag", "key"]

    def __str__(self):
        return self.key


class ValueParentMixin(Timestamped):
    value_type = models.ForeignKey(ValueType, on_delete=models.CASCADE)
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    object_id = models.PositiveIntegerField(
        help_text="The ID of the object in the table for that type of value (Float Value, Text Value, File Value)"
    )
    value_object = GenericForeignKey("content_type", "object_id")

    _value_models = ()

    class Meta:
        abstract = True

    @property
    def value(self):
        return self.value_object.value

    def clean(self):
        if self.content_type != self.value_type.content_type:
            raise ValidationError(
                _(
                    f"Value object's content type {self.content_type} doesn't "
                    "match value type's {self.value_type.content_type}"
                )
            )
        if self.content_type.model_class() not in self._value_models:
            raise ValidationError(
                _(f"Invalid model for solution: {self.content_type.model_class()}")
            )


class InputValue(ValueParentMixin):
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["input_element", "value_type"]

    def __str__(self):
        return f"{self.input_element}: {self.value_type}: {self.__str_value()}"

    def __str_value(self):
        # pylint: disable=no-member
        if isinstance(self.value, str) and len(self.value) > 100:
            return f"{self.value:.100}..."
        return str(self.value)

    def clean(self):
        if self.input_element.challenge.id != self.value_type.challenge.id:
            raise ValidationError(
                {"input_element": "Challenge must match value_type's"}
            )


class Evaluation(Logged):
    submission_run = models.ForeignKey(SubmissionRun, on_delete=models.CASCADE)
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)
    status = models.CharField(
        max_length=25, choices=Status.choices, default=Status.PENDING
    )

    class Meta:
        unique_together = ["submission_run", "input_element"]

    def mark_started(self, kwargs, file_kwargs):
        self.append(
            stdout="Started\n"
            f"Input element: {self.input_element}\n"
            f"kwargs: {kwargs}\n"
            f"file_kwargs: {file_kwargs}\n"
        )
        self.status = Status.RUNNING
        self.save()

    def __str__(self):
        return f"{self.submission_run}:, status {self.status}"

    def cleanup_local_outputs(self, output_file_keys):
        for prediction in self.prediction_set.filter(
            value_type__key__in=output_file_keys
        ):
            filecache.delete_local_cache(prediction.value)


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
        return self.key


class ScoreBase(Timestamped):
    score_type = models.ForeignKey(ScoreType, on_delete=models.CASCADE)
    value = models.FloatField()

    REQUIRED_LEVEL = None

    class Meta:
        abstract = True

    def clean(self):
        if self.score_type.level != self.REQUIRED_LEVEL:
            raise ValueError(
                f"Score Type {self.score_type} cannot be set on an {self.REQUIRED_LEVEL} score"
            )


class EvaluationScore(ScoreBase):
    evaluation = models.ForeignKey(
        Evaluation, on_delete=models.CASCADE, related_name="scores"
    )

    REQUIRED_LEVEL = ScoreType.Level.EVALUATION

    class Meta:
        unique_together = ["evaluation", "score_type"]

    def __str__(self):
        return f"{self.evaluation}:, {self.score_type}:{self.value}"


class SubmissionRunScore(ScoreBase):
    submission_run = models.ForeignKey(
        SubmissionRun, on_delete=models.CASCADE, related_name="scores"
    )

    REQUIRED_LEVEL = ScoreType.Level.SUBMISSION_RUN

    class Meta:
        unique_together = ["submission_run", "score_type"]

    def __str__(self):
        return f"{self.submission_run}:, {self.score_type}:{self.value}"


class Solution(ValueParentMixin):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)

    class Meta:
        abstract = True

    @classmethod
    def dicts_by_key(cls, instances):
        by_key = {}
        files_by_key = {}
        for instance in instances:
            if isinstance(instance.value_object, FileValue):
                files_by_key[instance.value_type.key] = filecache.ensure_local_copy(
                    instance.value
                )
            else:
                by_key[instance.value_type.key] = instance.value
        return by_key, files_by_key


class Prediction(Solution):
    evaluation = models.ForeignKey(Evaluation, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["evaluation", "value_type"]

    @classmethod
    def load_output(cls, challenge, evaluation, output_type, value):
        assert challenge is not None
        assert evaluation is not None
        assert output_type is not None

        prediction = cls(
            challenge=challenge,
            value_type=output_type,
            evaluation=evaluation,
        )
        output_type_model = output_type.content_type.model_class()
        value_object = output_type_model.from_string(
            value, challenge=challenge, evaluation=evaluation
        )

        value_object.save()
        prediction.value_object = value_object
        assert prediction.value_object is not None, "after save"

        return prediction

    def __str__(self):
        return f"{self.evaluation}::{self.value_type.key}::{self.content_type}"

    def clean(self):
        super().clean()
        self.challenge = self.evaluation.submission_run.submission.challenge


class AnswerKey(Solution):
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["input_element", "value_type"]

    def __str__(self):
        return f"{self.input_element}::{self.value_type.key}::{self.content_type}"

    def clean(self):
        super().clean()
        self.challenge = self.input_element.challenge


class GenericValue(models.Model):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    evaluation = models.ForeignKey(
        Evaluation, on_delete=models.CASCADE, null=True, blank=True
    )
    prediction = GenericRelation(Prediction)
    answer_key = GenericRelation(AnswerKey)
    input_value = GenericRelation(InputValue)

    class Meta:
        abstract = True

    @classmethod
    def from_string(cls, raw_value, *, challenge, evaluation=None, **_kwargs):
        cls_kwargs = {"challenge": challenge, "evaluation": evaluation}
        value_field = cls._meta.get_field("value")
        value = value_field.to_python(raw_value)
        return cls(value=value, **cls_kwargs)

    def __str__(self):
        # pylint: disable=no-member
        if isinstance(self.value, str) and len(self.value) > 100:
            return f"{self.value:.100}..."
        return str(self.value)


VALUE_PARENT_CLASSES = (Solution, InputValue)


def register_value_model(ValueModel):  # pylint: disable=invalid-name
    """
    Solution's value_object can only point to a
    value model registered with this decorator
    """
    if not hasattr(ValueModel, "value"):
        raise ValidationError(_(f"{ValueModel} must have a value attribute"))
    if not issubclass(ValueModel, GenericValue):
        raise ValidationError(_(f"{ValueModel} must extend GenericOutputValue"))
    for cls in VALUE_PARENT_CLASSES:
        if not issubclass(cls, ValueParentMixin):
            raise ValidationError(_(f"{cls} must extend ValueParentMixin"))
        cls._value_models = (  # pylint: disable=protected-access
            *cls._value_models,  # pylint: disable=protected-access
            ValueModel,
        )
    return ValueModel


@register_value_model
class TextValue(GenericValue):
    value = models.TextField(blank=True)


@register_value_model
class FloatValue(GenericValue):
    value = models.FloatField()


def _upload_location(instance, filename):
    challenge = instance.challenge
    evaluation = instance.evaluation
    if evaluation:
        parent_path = os.path.join("evaluations", f"{evaluation.id}")
    else:
        parent_path = os.path.join("challenges", f"{challenge.id}")
    return os.path.join("file_uploads", parent_path, filename)


@register_value_model
class FileValue(GenericValue):
    value = models.FileField(upload_to=_upload_location)

    @classmethod
    def from_string(
        cls, filepath, *, challenge, evaluation=None
    ):  # pylint:disable=arguments-renamed,arguments-differ
        cls_kwargs = {"challenge": challenge, "evaluation": evaluation}
        filename = os.path.basename(filepath)
        instance = cls(value=filename, **cls_kwargs)
        with open(filepath, "rb") as fp:
            instance.value.save(filename, File(fp))
        filecache.preserve_local_copy(instance.value, filepath)
        return instance
