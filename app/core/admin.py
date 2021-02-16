from django.contrib import admin
from django.contrib.admin import register

from . import models


class TimestampedAdmin(admin.ModelAdmin):
    readonly_fields = ("created_at", "updated_at")


# pylint: disable=no-self-use
# Sometimes we add methods to ModelAdmin that don't need to use self


@register(models.Challenge)
class ChallengeAdmin(TimestampedAdmin):
    list_display = ("name", "start_at", "end_at")
    date_hierarchy = "start_at"


@register(models.Container)
class ContainerAdmin(TimestampedAdmin):
    list_display = ("name", "user", "challenge", "created_at")
    date_hierarchy = "created_at"


@register(models.Submission)
class SubmissionAdmin(TimestampedAdmin):
    list_display = ("challenge", "user", "container", "created_at")


@register(models.SubmissionEvaluation)
class SubmissionEvaluationAdmin(TimestampedAdmin):
    list_display = (
        "submission",
        "user",
        "challenge",
        "exit_status",
    )
    date_hierarchy = "started_at"

    def user(self, instance):
        return instance.submission.user

    def challenge(self, instance):
        return instance.submission.challenge
