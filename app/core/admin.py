from django.contrib import admin
from django.contrib.admin import register
from django.contrib.admin.templatetags import admin_urls
from django.contrib.contenttypes.models import ContentType
from django.urls import reverse
from django.utils.html import format_html, format_html_join, mark_safe

from . import models

HREF_TEMPLATE = '<a href="{}">{}</a> {}'


def _admin_url(obj):
    return reverse(
        admin_urls.admin_urlname(
            obj._meta, "change"  #  pylint: disable=protected-access
        ),
        args=[obj.pk],
    )


def _admin_link(obj):
    return format_html(
        HREF_TEMPLATE, _admin_url(obj), obj, getattr(obj, "created_at", "")
    )


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
    readonly_fields = ["batch_status", "created_at", "updated_at"]

    def batch_status(self, instance):
        batch_group = instance.current_batch_group()
        if batch_group is None:
            return "No batches yet"
        return format_html(
            "{} batches, created {}",
            batch_group.inputbatch_set.count(),
            batch_group.created_at,
        )


@register(models.Container)
class ContainerAdmin(TimestampedAdmin):
    list_display = ("name", "user", "challenge", "created_at")
    list_filter = ("challenge",)
    date_hierarchy = "created_at"
    readonly_fields = ("submissions", *TimestampedAdmin.readonly_fields)

    def submissions(self, instance):
        return _admin_links(instance.submission_set.all())


@register(models.ContainerArg)
class ContainerArgAdmin(TimestampedAdmin):
    list_display = ("key", "container")


@register(models.ScoreMaker)
class ScoreMakerAdmin(TimestampedAdmin):
    list_display = ("challenge", "container")
    list_filter = ("challenge",)


@register(models.ScoreType)
class ScoreTypeAdmin(TimestampedAdmin):
    list_display = ("key", "level", "challenge")
    list_filter = ("challenge",)


@register(models.EvaluationScore)
class EvaluationScoreAdmin(TimestampedAdmin):
    list_display = ("submission_run", "input_element", "score_type", "value")
    list_filter = ("score_type__challenge",)


@register(models.SubmissionRunScore)
class SubmissionRunScoreAdmin(TimestampedAdmin):
    list_display = ("submission_run", "score_type", "value")
    list_filter = ("score_type__challenge",)


@register(models.Submission)
class SubmissionAdmin(TimestampedAdmin):
    list_display = ("challenge", "user", "container", "created_at")
    list_filter = ("challenge",)
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
    list_filter = (
        "submission__challenge",
        "status",
    )
    readonly_fields = (
        "submission",
        "evaluations",
        "batch_evaluations",
        "scores",
        *TimestampedAdmin.readonly_fields,
    )

    def user(self, instance):
        return instance.submission.user

    def challenge(self, instance):
        return instance.submission.challenge

    def evaluations(self, instance):
        return _admin_links(instance.evaluation_set.all())

    def batch_evaluations(self, instance):
        return _admin_links(instance.batchevaluation_set.all())

    def scores(self, instance):
        return _scores_list(instance)


@register(models.InputElement)
class InputElementAdmin(TimestampedAdmin):
    list_display = ("name", "challenge", "is_public")
    list_filter = ("challenge",)


@register(models.ValueType)
class ValueTypeAdmin(TimestampedAdmin):
    list_display = ("key", "challenge", "description", "content_type", "is_input_flag")
    list_filter = ("challenge",)


class ValueParentAdminMixin(TimestampedAdmin):
    readonly_fields = (
        "value_type_challenge",
        "input_element_challenge",
        "object_link",
        "value_object_challenge",
        *TimestampedAdmin.readonly_fields,
    )

    def input_element_challenge(self, instance):
        return instance.input_element.challenge

    def value_type_challenge(self, instance):
        return instance.value_type.challenge

    def value_object_challenge(self, instance):
        return instance.value_object.challenge

    def object_link(self, instance):
        return _admin_link(instance.value_object)

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == "content_type":
            allowed_models = (
                content_type.id
                for content_type in ContentType.objects.get_for_models(
                    *models.Solution._value_models  #  pylint: disable=protected-access
                ).values()
            )
            kwargs["queryset"] = ContentType.objects.filter(id__in=allowed_models)
        return super().formfield_for_foreignkey(db_field, request, **kwargs)


@register(models.InputValue)
class InputValueAdmin(ValueParentAdminMixin):
    list_display = (
        "pk",
        "input_element",
        "value_type",
        "content_type",
    )
    list_filter = ("input_element__challenge",)
    fields = (
        "download_file",
        ("value_type", "value_type_challenge"),
        ("content_type", "object_id"),
        ("object_link", "value_object_challenge"),
        ("input_element", "input_element_challenge"),
        ("created_at", "updated_at"),
    )
    readonly_fields = ("download_file", *ValueParentAdminMixin.readonly_fields)

    def download_file(self, instance):
        if instance.value_type.content_type.model != "filevalue":
            return "not a file"

        return format_html(
            HREF_TEMPLATE,
            reverse("download-input", args=[instance.pk]),
            instance.value_type.content_type,
            "",
        )

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        url_name = request.resolver_match.url_name
        if url_name == "core_inputvalue_change":
            input_value_id = int(request.resolver_match.kwargs["object_id"])
            input_value = models.InputValue.objects.get(pk=input_value_id)

            if db_field.name == "input_element":
                kwargs["queryset"] = (
                    input_value.value_type.challenge.inputelement_set.order_by("name")
                )

        return super().formfield_for_foreignkey(db_field, request, **kwargs)


@register(models.Evaluation)
class EvaluationAdmin(TimestampedAdmin):
    list_display = ("pk", "submission_run", "status")
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
        "submission_run",
        "input_element",
        "value_type",
        "content_type",
    )
    list_filter = ("challenge",)
    readonly_fields = ("challenge", "value", *TimestampedAdmin.readonly_fields)

    def value(self, instance):
        if instance.value_object:
            url = reverse(
                admin_urls.admin_urlname(
                    instance.value_object._meta,  #  pylint: disable=protected-access
                    "change",
                ),
                args=[instance.object_id],
            )
            print(url)
            return format_html('<a href="{}">{}</a>', url, instance.value_object)

        return "No value"


@register(models.AnswerKey)
class AnswerKeyAdmin(ValueParentAdminMixin):
    list_display = (
        "pk",
        "challenge",
        "input_element",
        "value_type",
        "content_type",
    )
    list_filter = ("challenge",)
    fields = (
        "challenge",
        ("value_type", "value_type_challenge"),
        ("content_type", "object_id"),
        ("object_link", "value_object_challenge"),
        ("input_element", "input_element_challenge"),
        ("created_at", "updated_at"),
    )

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        url_name = request.resolver_match.url_name
        if url_name == "core_answerkey_change":
            answer_key_id = int(request.resolver_match.kwargs["object_id"])
            answer_key = models.AnswerKey.objects.get(pk=answer_key_id)

            if db_field.name == "input_element":
                kwargs["queryset"] = answer_key.challenge.inputelement_set.order_by(
                    "name"
                )

        return super().formfield_for_foreignkey(db_field, request, **kwargs)


class GenericValueAdmin(admin.ModelAdmin):
    list_display = ("id", "__str__", "challenge", "evaluation")
    list_filter = ("challenge",)
    readonly_fields = ("prediction", "answer_key", "input_element")

    def input_element(self, instance):
        if instance.input_element.exists():
            url = reverse(
                "admin:core_inputelement_change", args=[instance.input_element.get().pk]
            )
            return format_html('<a href="{}">{}</a>', url, "Input Element")
        return "No input element"

    def prediction(self, instance):
        if instance.prediction.exists():
            url = reverse(
                "admin:core_prediction_change", args=[instance.prediction.get().pk]
            )
            return format_html('<a href="{}">{}</a>', url, "Prediction")
        return "No prediction"

    def answer_key(self, instance):
        if instance.answer_key.exists():
            url = reverse(
                "admin:core_answerkey_change", args=[instance.answer_key.get().pk]
            )
            return format_html('<a href="{}">{}</a>', url, "Answer Key")
        return "No answer key"


@register(models.TextValue)
class TextValueAdmin(GenericValueAdmin):
    pass


@register(models.FloatValue)
class FloatValueAdmin(GenericValueAdmin):
    pass


@register(models.FileValue)
class FileValueAdmin(GenericValueAdmin):
    pass


@register(models.InputBatchGroup)
class InputBatchGroupAdmin(TimestampedAdmin):
    list_filter = ("challenge",)
    list_display = ("id", "challenge", "created_at")


@register(models.InputBatch)
class InputBatchAdmin(TimestampedAdmin):
    list_display = ("id", "batch_group", "is_public")


@register(models.InputBatchMembership)
class InputBatchMembershipAdmin(TimestampedAdmin):
    list_display = ("id", "batch_group", "batch", "input_element")


@register(models.BatchFile)
class BatchFileAdmin(TimestampedAdmin):
    list_display = ("id", "batch", "value_type")
    readonly_fields = (
        "batch",
        "value_type",
        "download_batch",
        *TimestampedAdmin.readonly_fields,
    )

    def download_batch(self, instance):
        return format_html(
            HREF_TEMPLATE,
            reverse("download-batch", args=[instance.pk]),
            instance.data,
            "",
        )


@register(models.BatchEvaluation)
class BatchEvaluationAdmin(TimestampedAdmin):
    list_display = ("id", "input_batch", "status", "submission_run")
    readonly_fields = (
        "input_batch",
        "status",
        "submission_run",
        "log_stdout",
        "log_stderr",
        *TimestampedAdmin.readonly_fields,
    )
