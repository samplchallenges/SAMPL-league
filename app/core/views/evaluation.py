from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.views.generic.detail import DetailView

from core.models.run_related import BatchEvaluation

from ..models import Evaluation


class EvaluationUserTestMixin(LoginRequiredMixin, UserPassesTestMixin):
    def test_func(self):
        if self.request.user.is_superuser:
            return True
        evaluation = self.get_object()
        return evaluation.submission_run.submission.user == self.request.user


class BaseEvaluationDetail(EvaluationUserTestMixin, DetailView):
    context_object_name = "evaluation"

    def get_queryset(self):
        return (
            super()
            .get_queryset()
            .select_related("submission_run__submission__challenge")
        )

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["submission"] = context["evaluation"].submission_run.submission
        context["challenge"] = context["submission"].challenge
        return context


class EvaluationDetail(BaseEvaluationDetail):
    model = Evaluation


class BatchEvaluationDetail(BaseEvaluationDetail):
    model = BatchEvaluation


class BaseEvaluationLog(EvaluationUserTestMixin, DetailView):
    logtype = None
    template_name = "core/evaluation_log.html"

    def get(self, request, *args, **kwargs):
        self.logtype = kwargs.get("log")
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if self.logtype == "err":
            context["log"] = self.object.log_stderr
        elif self.logtype == "out":
            context["log"] = self.object.log_stdout
        else:
            raise ValueError(f"Invalid logtype: {self.logtype}")
        return context


class EvaluationLog(BaseEvaluationLog):
    model = Evaluation


class BatchEvaluationLog(BaseEvaluationLog):
    model = BatchEvaluation
