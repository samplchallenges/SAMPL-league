from django.db import models
import django.contrib.auth.models as auth_models
from django.utils import timezone


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
    # Or blobs? We can store pretty big files as blobs.
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

    def __str__(self):
        return str(self.name)


class Submission(Timestamped):
    user = models.ForeignKey(auth_models.User, on_delete=models.CASCADE)
    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    container = models.ForeignKey(Container, on_delete=models.CASCADE)
    doi = models.CharField(max_length=255)
    url = models.URLField()

    def __str__(self):
        return f"{self.user}: {self.challenge}"


class SubmissionField(Timestamped):
    name = models.CharField(max_length=255)
    description = models.TextField()


class SubmissionProperty(Timestamped):
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    field = models.ForeignKey(SubmissionField, on_delete=models.CASCADE)
    value = models.TextField()


class SubmissionEvaluation(Timestamped):
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    started_at = models.DateTimeField()
    ended_at = models.DateTimeField()
    exit_status = models.IntegerField()
    log_stdout = models.TextField(blank=True)
    log_stderr = models.TextField(blank=True)

    def __str__(self):
        return f"{self.submission}, exited {self.exit_status}"


class SubmissionResult(Timestamped):
    submission_evaluation = models.ForeignKey(
        SubmissionEvaluation, on_delete=models.CASCADE
    )
    # TBD: allow more than one result file per submission?
    datafile = models.FileField()
