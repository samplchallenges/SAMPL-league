from django.http import HttpResponseForbidden
from django.urls import reverse_lazy
from django.views.generic.edit import CreateView, DeleteView, UpdateView
from django.views.generic.detail import DetailView
from django.views.generic.list import ListView
from django.views.generic.detail import SingleObjectMixin

from ..forms import SubmissionForm
from ..models import Submission


class SubmissionDetail(DetailView):
    model = Submission


class SubmissionList(ListView):
    model = Submission


class SubmissionCreate(CreateView):
    model = Submission
    form_class = SubmissionForm


class SubmissionUpdate(UpdateView):
    model = Submission
    form_class = SubmissionForm


class SubmissionDelete(DeleteView):
    model = Submission
    success_url = reverse_lazy("submission-list")
