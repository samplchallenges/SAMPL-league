from django.contrib import admin
from django.contrib.admin import register
from django.contrib.admin.templatetags import admin_urls
from django.urls import reverse
from django.utils.html import format_html, format_html_join, mark_safe

from . import models

HREF_TEMPLATE = '<a href="{}">{}</a> {}t'


def _admin_url(obj):
    return reverse(
        admin_urls.admin_urlname(obj._meta, "change"),
        args=[obj.pk],
    )


def _admin_link(obj):
    return format_html(HREF_TEMPLATE, _admin_url(obj), obj, obj.created_at)


def _admin_links(objects):
    return format_html_join(
        mark_safe("<br/>"),
        HREF_TEMPLATE,
        [(_admin_url(obj), obj, obj.updated_at) for obj in objects],
    )


def _scores_list(instance):
    return format_html_join(
        mark_safe("<br/>"),
        "{}: {}",
        [(score.score_type.key, score.value) for score in instance.scores.all()],
    )


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
    readonly_fields = ("submissions", *TimestampedAdmin.readonly_fields)

    def submissions(self, instance):
        return _admin_links(instance.submission_set.all())


@register(models.ScoreMaker)
class ScoreMakerAdmin(TimestampedAdmin):
    list_display = ("challenge", "container")


@register(models.Submission)
class SubmissionAdmin(TimestampedAdmin):
    list_display = ("challenge", "user", "container", "created_at")

    readonly_fields = ("submission_runs", *TimestampedAdmin.readonly_fields)

    def submission_runs(self, instance):
        return _admin_links(instance.submissionrun_set.all())


@register(models.SubmissionRun)
class SubmissionRunAdmin(TimestampedAdmin):
    list_display = (
        "submission",
        "user",
        "challenge",
        "is_public",
        "status",
    )
    readonly_fields = (
        "submission",
        "evaluations",
        "scores",
        *TimestampedAdmin.readonly_fields,
    )

    def user(self, instance):
        return instance.submission.user

    def challenge(self, instance):
        return instance.submission.challenge

    def evaluations(self, instance):
        return _admin_links(instance.evaluation_set.all())

    def scores(self, instance):
        return _scores_list(instance)


@register(models.InputElement)
class InputElementAdmin(TimestampedAdmin):
    pass


@register(models.ValueType)
class ValueTypeAdmin(TimestampedAdmin):
    pass


@register(models.InputValue)
class InputValueAdmin(TimestampedAdmin):
    list_display = (
        "pk",
        "input_element",
        "value_type",
        "content_type",
    )


@register(models.Evaluation)
class EvaluationAdmin(TimestampedAdmin):
    list_display = ("pk", "submission_run")
    readonly_fields = (
        "submission_run",
        "predictions",
        "scores",
        *TimestampedAdmin.readonly_fields,
    )

    def submission_run(self, instance):
        return _admin_link(instance.submission_run)

    def predictions(self, instance):
        return _admin_links(instance.prediction_set.all())

    def scores(self, instance):
        return _scores_list(instance)


@register(models.Prediction)
class PredictionAdmin(TimestampedAdmin):
    list_display = (
        "pk",
        "challenge",
        "evaluation",
        "value_type",
        "content_type",
    )
    readonly_fields = ("challenge", "value", *TimestampedAdmin.readonly_fields)

    def value(self, instance):
        if instance.value_object:
            url = reverse(
                admin_urls.admin_urlname(instance.value_object._meta, "change"),
                args=[instance.object_id],
            )
            print(url)
            return format_html('<a href="{}">{}</a>', url, instance.value_object)

        return "No value"


@register(models.AnswerKey)
class AnswerKeyAdmin(TimestampedAdmin):
    list_display = (
        "pk",
        "challenge",
        "input_element",
        "value_type",
        "content_type",
    )
    readonly_fields = ("challenge", "value", *TimestampedAdmin.readonly_fields)

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
    readonly_fields = ("prediction", "answer_key")

    def prediction(self, instance):
        if instance.prediction:
            url = reverse(
                "admin:core_prediction_change", args=[instance.prediction.get().pk]
            )
            return format_html('<a href="{}">{}</a>', url, "Prediction")
        return "No prediction"

    def answer_key(self, instance):
        if instance.answer_key:
            url = reverse(
                "admin:core_answerkey_change", args=[instance.answer_key.get().pk]
            )
            return format_html('<a href="{}">{}</a>', url, "Answer Key")
        return "No answer key"


@register(models.TextValue)
class TextValueAdmin(GenericOutputValueAdmin):
    pass


@register(models.FloatValue)
class FloatValueAdmin(GenericOutputValueAdmin):
    pass


@register(models.BlobValue)
class BlobValueAdmin(GenericOutputValueAdmin):
    pass
