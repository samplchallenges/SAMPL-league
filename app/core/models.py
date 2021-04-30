import django.contrib.auth.models as auth_models
from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ValidationError
from django.db import models
from django.urls import reverse
from django.utils import timezone
from django.utils.translation import gettext_lazy as _

from . import configurator


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

    def __str__(self):
        return str(self.name)


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
    class _Status(models.TextChoices):
        FAILURE = "FAILURE"
        SUCCESS = "SUCCESS"
        PENDING = "PENDING"

    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    digest = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)
    status = models.CharField(max_length=25, choices=_Status.choices)

    def __str__(self):
        return f"{self.submission}:{self.digest}, status {self.status}"


class InputElement(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    name = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)

    class Meta:
        unique_together = ["challenge", "name"]

    def __str__(self):
        return f"{self.name}, is public? {self.is_public}"


class InputType(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    key = models.CharField(max_length=255)
    description = models.TextField()

    class Meta:
        unique_together = ["challenge", "key"]

    def __str__(self):
        return self.key


class InputValue(Timestamped):
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)
    input_type = models.ForeignKey(InputType, on_delete=models.CASCADE)
    value = models.TextField()
    # TODO: use https://docs.djangoproject.com/en/3.2/ref/contrib/contenttypes/

    class Meta:
        unique_together = ["input_element", "input_type"]

    def __str__(self):
        return f"{self.input_element}: {self.input_type}: {self.__str_value()}"

    def __str_value(self):
        # pylint: disable=no-member
        if isinstance(self.value, str) and len(self.value) > 100:
            return f"{self.value:.100}..."
        return str(self.value)


class Evaluation(Timestamped):
    submission_run = models.ForeignKey(SubmissionRun, on_delete=models.CASCADE)
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)
    exit_status = models.IntegerField()
    log_stdout = models.TextField(blank=True)
    log_stderr = models.TextField(blank=True)

    class Meta:
        unique_together = ["submission_run", "input_element"]

    def __str__(self):
        return f"{self.submission_run}:, exited {self.exit_status}"


class Solution(Timestamped):
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    key = models.CharField(max_length=255)
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    object_id = models.PositiveIntegerField()
    value_object = GenericForeignKey("content_type", "object_id")

    _value_models = ()

    class Meta:
        abstract = True

    @classmethod
    def register_value_model(cls, ValueModel):
        """
        Solution's value_object can only point to a
        value model registered with this decorator
        """
        if not hasattr(ValueModel, "value"):
            raise ValidationError(_(f"{ValueModel} must have a value attribute"))
        if not issubclass(ValueModel, GenericOutputValue):
            raise ValidationError(_(f"{ValueModel} must extend GenericOutputValue"))

        cls._value_models = (*cls._value_models, ValueModel)
        return ValueModel

    def clean(self):
        if self.content_type.model_class() not in self._value_models:
            raise ValidationError(
                _(f"Invalid model for solution: {self.content_type.model_class()}")
            )


class Prediction(Solution):
    evaluation = models.ForeignKey(Evaluation, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["evaluation", "key"]

    def __str__(self):
        return f"{self.evaluation}::{self.key}::{self.content_type}"

    def clean(self):
        super().clean()
        self.challenge = self.evaluation.submission_run.submission.challenge


class AnswerKey(Solution):
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["input_element", "key"]

    def __str__(self):
        return f"{self.input_element}::{self.key}::{self.content_type}"

    def clean(self):
        super().clean()
        self.challenge = self.input_element.challenge


class GenericOutputValue(models.Model):
    prediction = GenericRelation(Prediction)

    class Meta:
        abstract = True

    def __str__(self):
        # pylint: disable=no-member
        if isinstance(self.value, str) and len(self.value) > 100:
            return f"{self.value:.100}..."
        return str(self.value)


@Solution.register_value_model
class TextValue(GenericOutputValue):
    value = models.TextField(blank=True)


@Solution.register_value_model
class FloatValue(GenericOutputValue):
    value = models.FloatField()


@Solution.register_value_model
class BlobValue(GenericOutputValue):
    value = models.BinaryField()
