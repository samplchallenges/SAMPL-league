from django.db import models
from django.utils import timezone


class ContainerType(models.TextChoices):
    DOCKER = "docker"
    SINGULARITY_LOCAL = "singularity_local"
    # TODO: Singularity remote server


class Timestamped(models.Model):
    created_at = models.DateTimeField(default=timezone.now)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        abstract = True
        get_latest_by = ("created_at",)
