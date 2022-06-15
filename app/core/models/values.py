import os.path

from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ValidationError
from django.core.files import File
from django.db import models
from django.utils.translation import gettext_lazy as _

from .. import filecache
from .admin_managed import Challenge, InputElement, ValueType
from .infra_models import Timestamped


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

        if self.input_element.is_parent and not self.value_type.on_parent_flag:
            raise ValidationError(
                "Cannot set a non-parent value type on a parent input element"
            )

        if not self.input_element.is_parent and self.value_type.on_parent_flag:
            raise ValidationError(
                "Cannot set a parent value type on non-parent input element"
            )


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
        "Evaluation", on_delete=models.CASCADE, null=True, blank=True
    )
    batch_evaluation = models.ForeignKey(
        "BatchEvaluation", on_delete=models.CASCADE, null=True, blank=True
    )
    prediction = GenericRelation("Prediction")
    batch_prediction = GenericRelation("BatchPrediction")

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
