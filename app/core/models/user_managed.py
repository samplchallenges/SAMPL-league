import os.path

import django.contrib.auth.models as auth_models
from django.core.exceptions import ValidationError
from django.db import models
from django.urls import reverse

from .. import configurator
from .admin_managed import Challenge, Container
from .infra_models import Timestamped


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

    def create_run(self, *, is_public, status):
        # TODO: move this method into referee's tasks.py?
        # status = Status.PENDING_REMOTE if remote else Status.PENDING
        digest = self.container.digest
        if digest is None:
            digest = "nodigest"
        submission_run = self.submissionrun_set.create(
            digest=digest,
            is_public=is_public,
            status=status,
        )

        batch_group = self.challenge.current_batch_group()
        if batch_group is None:
            for element_id in self.challenge.inputelement_set.filter(
                is_public=is_public
            ).values_list("id", flat=True):
                submission_run.evaluation_set.create(input_element_id=element_id)
        else:
            for input_batch_id in batch_group.inputbatch_set.filter(
                is_public=is_public
            ).values_list("id", flat=True):
                submission_run.batchevaluation_set.create(input_batch_id=input_batch_id)

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
