from django.utils.translation import gettext_lazy as _
from django.db import models
import django.contrib.auth.models as auth_models
from django.utils import timezone
from django.urls import reverse

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
    # TODO: How to store scoring container reference in challenge?

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
                if isinstance(field, models.fields.Field):
                    if not getattr(self, field.attname):
                        self.draft_mode = True
                        break


class SubmissionEvaluation(Timestamped):
    class _DataPrivacyLevel(models.TextChoices):
        PUBLIC = "PUBLIC"
        PRIVATE = "PRIVATE"

    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    digest = models.CharField(max_length=255)
    data_privacy_level = models.CharField(
        max_length=25,
        choices=_DataPrivacyLevel.choices,
        default=_DataPrivacyLevel.PRIVATE,
    )
    started_at = models.DateTimeField()
    ended_at = models.DateTimeField()
    exit_status = models.IntegerField()
    log_stdout = models.TextField(blank=True)
    log_stderr = models.TextField(blank=True)

    def __str__(self):
        return f"{self.submission}:{self.digest}, exited {self.exit_status}"


class SubmissionResult(Timestamped):
    submission_evaluation = models.ForeignKey(
        SubmissionEvaluation, on_delete=models.CASCADE
    )
    # TBD: allow more than one result file per submission?
    datafile = models.FileField()
