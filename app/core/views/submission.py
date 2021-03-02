from django.http import HttpResponseForbidden, HttpResponseBadRequest
from django.shortcuts import redirect, render
from django.urls import reverse_lazy
from django.views.generic.edit import CreateView, DeleteView, UpdateView
from django.views.generic.detail import DetailView
from django.views.generic.list import ListView
from django.views.generic.detail import SingleObjectMixin

from ..forms import SubmissionForm, ContainerForm
from ..models import Submission


class SubmissionDetail(DetailView):
    model = Submission


class SubmissionList(ListView):
    model = Submission


# class SubmissionCreate(CreateView):
#    model = Submission
#    form_class = SubmissionForm


# class SubmissionUpdate(UpdateView):
#    model = Submission
#    form_class = SubmissionForm


class SubmissionDelete(DeleteView):
    model = Submission
    success_url = reverse_lazy("submission-list")


def edit_submission_view(request):
    if request.method == "POST":
        # TODO: transaction?
        container_form = ContainerForm(request.POST)
        submission_form = SubmissionForm(request.POST)
        if container_form.is_valid():
            container = container_form.save(commit=False)
            container.user = request.user
            container.save()
            submission_form.container = container
            if submission_form.is_valid():
                submission = submission_form.save(commit=False)
                submission.user = request.user
                submission.save()
                return redirect("submission-detail", pk=submission.pk)
    elif request.method == "GET":
        container_form = ContainerForm()
        submission_form = SubmissionForm()
    else:
        return HttpResponseBadRequest()
    context = {"container_form": container_form, "submission_form": submission_form}
    return render(request, "core/submission_form.html", context)
