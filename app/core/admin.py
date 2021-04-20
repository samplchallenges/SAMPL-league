from django.contrib import admin
from django.contrib.admin import register
from django.contrib.admin.templatetags import admin_urls
from django.urls import reverse
from django.utils.html import format_html

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


@register(models.ScoreMaker)
class ScoreMakerAdmin(TimestampedAdmin):
    list_display = ("challenge", "container")


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


@register(models.InputElement)
class InputElementAdmin(TimestampedAdmin):
    pass


@register(models.InputType)
class InputTypeAdmin(TimestampedAdmin):
    pass


@register(models.InputValue)
class InputValueAdmin(TimestampedAdmin):
    pass


@register(models.Prediction)
class PredictionAdmin(TimestampedAdmin):
    list_display = (
        "pk",
        "submission_evaluation",
        "input_element",
        "key",
        "content_type",
    )
    readonly_fields = ("value",)

    def value(self, instance):
        if instance.value_object:
            url = reverse(
                admin_urls.admin_urlname(instance.value_object._meta, "change"),
                args=[instance.object_id],
            )
            print(url)
            return format_html('<a href="{}">{}</a>', url, instance.value_object)

        return "No value"


class GenericOutputValueAdmin(admin.ModelAdmin):
    readonly_fields = ("prediction",)

    def prediction(self, instance):
        if instance.prediction:
            url = reverse(
                "admin:core_prediction_change", args=[instance.prediction.get().pk]
            )
            return format_html('<a href="{}">{}</a>', url, "Prediction")

        return "No prediction"


@register(models.TextValue)
class TextValueAdmin(GenericOutputValueAdmin):
    pass


@register(models.FloatValue)
class FloatValueAdmin(GenericOutputValueAdmin):
    pass


@register(models.BlobValue)
class BlobValueAdmin(GenericOutputValueAdmin):
    pass
