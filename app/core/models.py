import logging
import os.path

import django.contrib.auth.models as auth_models
from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ValidationError
from django.core.files import File
from django.db import models
from django.urls import reverse
from django.utils import timezone
from django.utils.translation import gettext_lazy as _

from . import configurator
from . import filecache


logger = logging.getLogger(__name__)


class Status(models.TextChoices):
    FAILURE = "FAILURE"
    SUCCESS = "SUCCESS"
    PENDING = "PENDING"


class Timestamped(models.Model):
    created_at = models.DateTimeField(default=timezone.now)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        abstract = True
        get_latest_by = ("created_at",)


class Challenge(Timestamped):
    name = models.CharField(max_length=255)
    start_at = models.DateTimeField()
    end_at = models.DateTimeField()
    repo_url = models.URLField()
    # Data stored in S3 - privileges on S3 buckets will help
    # prevent leakage of secret data
    sample_data_url = models.URLField()
    sample_score_reference_url = models.URLField()
    secret_data_url = models.URLField()
    secret_score_reference_url = models.URLField()
    execution_options_json = models.JSONField()

    __output_types_dict = None

    def __str__(self):
        return str(self.name)

    def __load_output_types(self):
        output_types = self.valuetype_set.filter(is_input_flag=False)
        self.__output_types_dict = {
            output_type.key: output_type for output_type in output_types.all()
        }
        file_content_type = ContentType.objects.get_for_model(FileValue)
        self.__output_file_keys = set(
            key
            for key, output_type in self.__output_types_dict.items()
            if output_type.content_type == file_content_type
        )

    def output_file_keys(self):
        if self.__output_types_dict is None:
            self.__load_output_types()
        return self.__output_file_keys

    def output_type(self, key):
        if self.__output_types_dict is None:
            self.__load_output_types()
        return self.__output_types_dict.get(key)


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

    def __str__(self):
        return f"{self.user}: {self.challenge}: {self.name}"

    def get_absolute_url(self):
        return reverse("submission", kwargs={"pk": self.pk})

    def clean(self):
        super().clean()
        if not self.draft_mode:
            # All fields
            for field in self._meta.get_fields(include_parents=False):
                if isinstance(field, models.fields.Field) and not isinstance(
                    field, models.BooleanField
                ):
                    if not getattr(self, field.attname):
                        self.draft_mode = True
                        break

    def clone(self):
        """
        Return a new submission using this as a template.
        Set draft mode and clear 'Name' field.
        """
        self.id = None
        self.pk = None
        self._state.adding = True
        return self


class SubmissionRun(Timestamped):
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    digest = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)
    status = models.CharField(max_length=25, choices=Status.choices)

    def __str__(self):
        return f"{self.submission}:{self.digest}, status {self.status}"


class InputElement(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    name = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)

    class Meta:
        unique_together = ["challenge", "name"]

    def all_values(self):
        """
        Returns a pair of key: value dicts, where the first dict is regular values
        and the second is file values
        """
        file_values = {}
        values = {}
        for input_value in self.inputvalue_set.select_related(
            "value_type", "value_type__content_type"
        ).all():
            value_type = input_value.value_type
            key = value_type.key
            content_type = value_type.content_type
            if content_type.model_class() == FileValue:
                file_values[key] = filecache.ensure_local_copy(input_value.value)
            else:
                values[key] = input_value.value
        return values, file_values

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
    object_id = models.PositiveIntegerField()
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


class Evaluation(Timestamped):
    submission_run = models.ForeignKey(SubmissionRun, on_delete=models.CASCADE)
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)
    status = models.CharField(
        max_length=25, choices=Status.choices, default=Status.PENDING
    )
    log_stdout = models.TextField(blank=True)
    log_stderr = models.TextField(blank=True)

    class Meta:
        unique_together = ["submission_run", "input_element"]

    def __str__(self):
        return f"{self.submission_run}:, status {self.status}"


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
            if isinstance(instance.value_type, FileValue):
                files_by_key[instance.value_type.key] = instance.value.path
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
        logger.debug(f"{prediction.value_object.__dict__}")

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


def register_value_model(ValueModel):
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
        cls._value_models = (*cls._value_models, ValueModel)
    return ValueModel


@register_value_model
class TextValue(GenericValue):
    value = models.TextField(blank=True)


@register_value_model
class FloatValue(GenericValue):
    value = models.FloatField()


def _upload_location(instance, filename):
    now = timezone.now()
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
    def from_string(cls, filepath, *, challenge, evaluation=None):
        cls_kwargs = {"challenge": challenge, "evaluation": evaluation}
        filename = os.path.basename(filepath)
        instance = cls(value=filename, **cls_kwargs)
        with open(filepath) as fp:
            instance.value.save(filename, File(fp))
        return instance
