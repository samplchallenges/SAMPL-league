import os.path
import time
from collections import defaultdict, namedtuple
from functools import partial

from django.db import models, transaction
from django.db.models.functions import Concat
from ever_given.utils import LogHandlerBase

from core import values_helper

from .. import batching, filecache
from .admin_managed import InputElement, ScoreType
from .batch_related import BatchFileBase, InputBatch
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

    def update_status(self, status):
        self.status = status
        self.save(update_fields=["status"])


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

    def mark_for_cancel(self, remote=False):
        with transaction.atomic():
            Evaluation.objects.select_for_update().select_related(
                "submission_run"
            ).filter(submission_run_id=self.id).all()
            submission_run = SubmissionRun.objects.select_for_update().get(pk=self.id)
            # TODO: handle batches in cancel
            if remote and submission_run.status == Status.PENDING_REMOTE:
                self.update_status(status=Status.CANCELLED)
                self.evaluation_set.update(status=Status.CANCELLED)
            else:
                self.update_status(status=Status.CANCEL_PENDING)
                self.evaluation_set.filter(status=Status.PENDING).update(
                    status=Status.CANCELLED
                )
                self.evaluation_set.filter(status=Status.RUNNING).update(
                    status=Status.CANCEL_PENDING
                )

    def completion(self):
        if self.submission.challenge.max_batch_size > 0:
            batch_group = self.submission.challenge.current_batch_group()
            completed = self.batchevaluation_set.filter(
                status__in=(Status.SUCCESS, Status.FAILURE)
            ).count()
            num_runs = batch_group.inputbatch_set.filter(
                is_public=self.is_public
            ).count()
        else:
            completed = self.evaluation_set.filter(
                status__in=(Status.SUCCESS, Status.FAILURE)
            ).count()
            num_runs = self.submission.challenge.inputelement_set.filter(
                is_public=self.is_public, is_parent=False
            ).count()
        completed_frac = completed / num_runs if num_runs else 0
        return Completion(completed, num_runs, completed_frac)


class SubmissionRunPair(Timestamped):
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    public_run = models.ForeignKey(
        SubmissionRun, on_delete=models.CASCADE, related_name="%(class)s_public_run"
    )
    private_run = models.ForeignKey(
        SubmissionRun, on_delete=models.CASCADE, related_name="%(class)s_private_run"
    )

    class Meta:
        unique_together = ["public_run", "private_run"]

    def __str__(self):
        return (
            f"{self.submission}: public {self.public_run}, private {self.private_run}"
        )


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

    def clear_old_predictions(self, output_type):
        raise Exception("Must implement in subclass")

    def clear_old_scores(self):
        self.scores.delete()

    def mark_started(self):
        self.append(stdout="Started\n")
        self.status = Status.RUNNING
        self.save(update_fields=["status"])

    def __str__(self):
        return f"run: {self.submission_run}, local status {self.status}"


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

    def clear_old_predictions(self, output_type):
        matching_prediction = Prediction.objects.filter(
            submission_run_id=self.submission_run_id,
            input_element_id=self.input_element_id,
            value_type=output_type,
        )
        if matching_prediction.exists():
            self.append(stderr="Duplicate prediction entry, overwriting old entry")
            matching_prediction.delete()

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
    def _build_single(
        cls, submission_run, value_type, input_element_id, value, **evaluation_attrs
    ):
        challenge = value_type.challenge
        prediction = cls(
            challenge=challenge,
            submission_run=submission_run,
            input_element_id=input_element_id,
            value_type=value_type,
        )
        output_type_model = value_type.content_type.model_class()
        if len(evaluation_attrs) != 1 or set(evaluation_attrs.keys()).isdisjoint(
            {"evaluation", "batch_evaluation"}
        ):
            raise ValueError(
                f"Must pass either evaluation or batch_evaluation, not {evaluation_attrs}"
            )
        value_object = output_type_model.from_string(
            value,
            challenge=challenge,
            input_element_id=input_element_id,
            **evaluation_attrs,
        )

        value_object.save()
        prediction.value_object = value_object
        if prediction.value_object is None:
            raise ValueError("after save, value_object must not be none")

        prediction.save()

    @classmethod
    def load_evaluation_output(cls, challenge, evaluation, output_type, value):
        if challenge is None:
            raise ValueError("Challenge is not set")
        if evaluation is None:
            raise ValueError("Evaluation is not set")
        if output_type is None:
            raise ValueError("Output type is not set")
        cls._build_single(
            submission_run=evaluation.submission_run,
            value_type=output_type,
            input_element_id=evaluation.input_element_id,
            value=value,
            evaluation=evaluation,
        )

    @classmethod
    def load_batch_output(cls, challenge, batch_evaluation, output_type, value):
        if challenge is None:
            raise ValueError("Challenge is not set")
        if batch_evaluation is None:
            raise ValueError("Batch Evaluation is not set")
        if output_type is None:
            raise ValueError("Output type is not set")
        if output_type.on_parent_flag:
            parent_elem_id = batch_evaluation.input_batch.parent_input_element_id
            cls._build_single(
                submission_run=batch_evaluation.submission_run,
                value_type=output_type,
                input_element_id=parent_elem_id,
                value=value,
                batch_evaluation=batch_evaluation,
            )
        else:
            batch_method = output_type.batch_method
            batch_cls = batching.BATCHERS[batch_method]
            for input_element_id, unbatched_value in batch_cls.invert(value):
                cls._build_single(
                    submission_run=batch_evaluation.submission_run,
                    value_type=output_type,
                    input_element_id=input_element_id,
                    value=unbatched_value,
                    batch_evaluation=batch_evaluation,
                )

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

    def __str__(self):
        visibility = "Public" if self.submission_run.is_public else "Private"
        num_elements = self.input_batch.inputbatchmembership_set.count()
        return f"{visibility} Batch #{self.input_batch.id} ({num_elements} elements)"

    def scores_dicts(self):
        by_elem = defaultdict(dict)
        for score in self.scores:
            by_elem[score.input_element_id][score.score_type.key] = score.value
        return by_elem.values()

    @property
    def input_object(self):
        return self.input_batch

    def predictions(self):
        return Prediction.objects.filter(
            submission_run_id=self.submission_run_id,
            input_element_id__in=[el.id for el in self.input_batch.elements()],
        )

    def clear_old_predictions(self, output_type):
        matching_objects = Prediction.objects.filter(
            submission_run_id=self.submission_run_id,
            input_element_id__in=[el.id for el in self.input_batch.elements()],
            value_type=output_type,
        )
        if matching_objects.exists():
            self.append(stderr="Duplicate prediction entry, overwriting old entry")
            matching_objects.delete()

    def cleanup_local_outputs(self, output_file_keys):
        for prediction in Prediction.objects.filter(
            submission_run=self.submission_run,
            input_element__in=self.input_batch.elements(),
            value_type__key__in=output_file_keys,
        ):
            filecache.delete_local_cache(prediction.value)

    def batchup(self):
        create_file_func = partial(
            BatchEvaluationFile.from_local, batch_evaluation=self
        )

        with transaction.atomic():
            value_types = self.input_batch.batch_group.challenge.valuetype_set.filter(
                is_input_flag=False, on_parent_flag=False
            )
            elements = list(self.input_batch.elements())

            batching.batchup_elements(
                value_types, elements, create_file_func, values_helper.predicted_values
            )


def batch_evaluation_upload_location(instance, filename):
    batch = instance.batch_evaluation.input_batch
    batch_group = batch.batch_group
    challenge = batch_group.challenge
    challenge_path = os.path.join("batch_evaluations", f"{challenge.id}")
    group_path = os.path.join("batch_group", f"{batch_group.id}")
    srun_path = os.path.join(
        "submission_run", f"{instance.batch_evaluation.submission_run.id}"
    )
    batch_path = os.path.join("b", f"{batch.id}")
    return os.path.join(
        challenge_path,
        group_path,
        srun_path,
        batch_path,
        instance.value_type.key,
        filename,
    )


class BatchEvaluationFile(BatchFileBase):
    batch_evaluation = models.ForeignKey(BatchEvaluation, on_delete=models.CASCADE)
    data = models.FileField(upload_to=batch_evaluation_upload_location)

    def __str__(self):
        return f"{self.batch_evaluation} {self.value_type} {self.data}"
