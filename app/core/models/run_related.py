import time
from collections import namedtuple

from django.db import models, transaction
from django.db.models.functions import Concat
from ever_given.utils import LogHandlerBase

from .. import filecache
from ..batching import BATCHERS
from .admin_managed import InputElement, ScoreType
from .batch_related import InputBatch
from .infra_models import Timestamped
from .user_managed import Submission
from .values import Solution

Completion = namedtuple("Completion", ["completed", "not_completed", "completed_frac"])


class Status(models.TextChoices):
    FAILURE = "FAILURE"
    SUCCESS = "SUCCESS"
    PENDING = "PENDING"
    PENDING_REMOTE = "PENDING_REMOTE"
    RUNNING = "RUNNING"
    CANCELLED = "CANCELLED"
    CANCEL_PENDING = "CANCEL_PENDING"


class StatusMixin:
    def is_finished(self):
        return self.status in (Status.SUCCESS, Status.FAILURE, Status.CANCELLED)


class Logged(Timestamped, StatusMixin):
    log_stdout = models.TextField(blank=True)
    log_stderr = models.TextField(blank=True)

    log_handler = None

    class Meta:
        abstract = True

    class LogHandler(LogHandlerBase):
        def __init__(self, instance):
            self.instance_id = instance.pk
            self.cls = type(instance)

        @property
        def _update(self):
            return self.cls.objects.filter(pk=self.instance_id).update

        def handle_stdout(self, log):
            self._update(
                log_stdout=Concat(
                    models.F("log_stdout"), models.Value(_timestamped_log(log))
                )
            )

        def handle_stderr(self, log):
            self._update(
                log_stderr=Concat(
                    models.F("log_stderr"), models.Value(_timestamped_log(log))
                )
            )

    def append(self, stdout=None, stderr=None):
        update_fields = []
        if stdout is not None:
            self.log_stdout = Concat(models.F("log_stdout"), models.Value(stdout))
            update_fields.append("log_stdout")
        if stderr is not None:
            self.log_stderr = Concat(models.F("log_stderr"), models.Value(stderr))
            update_fields.append("log_stderr")

        if update_fields:
            self.save(update_fields=update_fields)


class SubmissionRun(Logged):
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    digest = models.CharField(max_length=255)
    is_public = models.BooleanField(default=False)
    status = models.CharField(max_length=25, choices=Status.choices)

    def __str__(self):
        return f"{self.submission}:{self.digest}, status {self.status}"

    def check_cancel_requested(self):
        status = SubmissionRun.objects.values_list("status", flat=True).get(pk=self.id)
        return status == Status.CANCEL_PENDING

    def mark_for_cancel(self):
        with transaction.atomic():
            SubmissionRun.objects.filter(pk=self.id).update(
                status=Status.CANCEL_PENDING
            )
            Evaluation.objects.filter(
                submission_run_id=self.id, status=Status.PENDING
            ).update(status=Status.CANCELLED)
            Evaluation.objects.filter(
                submission_run_id=self.id, status=Status.RUNNING
            ).update(status=Status.CANCEL_PENDING)

    def completion(self):
        completed = self.evaluation_set.filter(
            status__in=(Status.SUCCESS, Status.FAILURE)
        ).count()
        num_element_ids = self.submission.challenge.inputelement_set.filter(
            is_public=self.is_public
        ).count()
        completed_frac = completed / num_element_ids if num_element_ids else 0
        return Completion(completed, num_element_ids, completed_frac)


class BaseEvaluation(Logged):
    submission_run = models.ForeignKey(SubmissionRun, on_delete=models.CASCADE)
    status = models.CharField(
        max_length=25, choices=Status.choices, default=Status.PENDING
    )

    class Meta:
        abstract = True

    @property
    def input_object(self):
        raise Exception("Must implement in subclasses")

    def mark_started(self, kwargs, file_kwargs):
        self.append(
            stdout="Started\n"
            f"Input object: {self.input_object}\n"
            f"kwargs: {kwargs}\n"
            f"file_kwargs: {file_kwargs}\n"
        )
        self.status = Status.RUNNING
        self.save()

    def __str__(self):
        return f"{self.submission_run}:, status {self.status}"


class Evaluation(BaseEvaluation):
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["submission_run", "input_element"]

    @property
    def scores(self):
        return EvaluationScore.objects.filter(
            submission_run_id=self.submission_run_id,
            input_element_id=self.input_element_id,
        )

    @property
    def input_object(self):
        return self.input_element

    def cleanup_local_outputs(self, output_file_keys):
        # TODO: copied between this and BatchEvaluation with minor changes
        for prediction in Prediction.objects.filter(
            submission_run=self.submission_run,
            input_element=self.input_element,
            value_type__key__in=output_file_keys,
        ):
            filecache.delete_local_cache(prediction.value)


class ScoreBase(Timestamped):
    submission_run = models.ForeignKey(SubmissionRun, on_delete=models.CASCADE)
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


class SubmissionRunScore(ScoreBase):
    REQUIRED_LEVEL = ScoreType.Level.SUBMISSION_RUN

    class Meta:
        unique_together = ["submission_run", "score_type"]

    def __str__(self):
        return f"{self.submission_run}:, {self.score_type}:{self.value}"


class EvaluationScore(ScoreBase):
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)

    REQUIRED_LEVEL = ScoreType.Level.EVALUATION

    class Meta:
        unique_together = ["submission_run", "input_element", "score_type"]

    def __str__(self):
        return f"{self.submission_run}::{self.input_element}:, {self.score_type}:{self.value}"


def _timestamped_log(log):
    return " ".join([time.strftime("[%c %Z]", time.gmtime()), log])


class Prediction(Solution):
    submission_run = models.ForeignKey(
        SubmissionRun, on_delete=models.CASCADE, related_name="predictions"
    )
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["submission_run", "input_element", "value_type"]

    @classmethod
    def load_evaluation_output(cls, challenge, evaluation, output_type, value):
        if challenge is None:
            raise ValueError("Challenge is not set")
        if evaluation is None:
            raise ValueError("Evaluation is not set")
        if output_type is None:
            raise ValueError("Output type is not set")

        prediction = cls(
            challenge=challenge,
            submission_run=evaluation.submission_run,
            input_element=evaluation.input_element,
            value_type=output_type,
        )
        output_type_model = output_type.content_type.model_class()
        value_object = output_type_model.from_string(
            value, challenge=challenge, input_element=evaluation.input_element
        )

        value_object.save()
        prediction.value_object = value_object
        if prediction.value_object is None:
            raise ValueError("after save, value_object must not be none")

        prediction.save()

    @classmethod
    def load_batch_output(cls, challenge, batch_evaluation, output_type, value):
        if challenge is None:
            raise ValueError("Challenge is not set")
        if batch_evaluation is None:
            raise ValueError("Batch Evaluation is not set")
        if output_type is None:
            raise ValueError("Output type is not set")

        batch_method = output_type.batch_method
        batch_cls = BATCHERS[batch_method]
        for input_element_id, unbatched_value in batch_cls.invert(value):
            prediction = cls(
                challenge=challenge,
                value_type=output_type,
                input_element_id=input_element_id,
                batch_evaluation=batch_evaluation,
            )
            output_type_model = output_type.content_type.model_class()
            value_object = output_type_model.from_string(
                unbatched_value, challenge=challenge, input_element_id=input_element_id
            )
            value_object.save()
            prediction.value_object = value_object
            if prediction.value_object is None:
                raise ValueError("after save, value_object must not be none")

            prediction.save()

    def __str__(self):
        return f"{self.input_element}::{self.value_type.key}::{self.content_type}"

    def clean(self):
        super().clean()
        self.challenge = self.submission_run.submission.challenge


class BatchEvaluation(BaseEvaluation):
    input_batch = models.ForeignKey(InputBatch, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["submission_run", "input_batch"]

    @property
    def scores(self):
        return EvaluationScore.objects.filter(
            submission_run_id=self.submission_run_id,
            input_element__in=self.input_batch.elements(),
        )

    @property
    def input_object(self):
        return self.input_batch

    def cleanup_local_outputs(self, output_file_keys):
        for prediction in Prediction.objects.filter(
            submission_run=self.submission_run,
            input_element__in=self.input_batch.elements(),
            value_type__key__in=output_file_keys,
        ):
            filecache.delete_local_cache(prediction.value)
