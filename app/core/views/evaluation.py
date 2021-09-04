from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.views.generic.detail import DetailView

from ..models import Evaluation


class EvaluationUserTestMixin(LoginRequiredMixin, UserPassesTestMixin):
    def test_func(self):
        if self.request.user.is_superuser:
            return True
        evaluation = self.get_object()
        return evaluation.submission_run.submission.user == self.request.user


class EvaluationDetail(EvaluationUserTestMixin, DetailView):
    model = Evaluation


class EvaluationLog(EvaluationUserTestMixin, DetailView):
    model = Evaluation
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
            raise ValueError("Invalid logtype: %s", self.logtype)
        return context
