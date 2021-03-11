from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.auth.mixins import UserPassesTestMixin
from django.http import HttpResponseForbidden, HttpResponseBadRequest
from django.shortcuts import redirect, render
from django.urls import reverse_lazy
from django.views.generic.edit import CreateView, DeleteView, UpdateView
from django.views.generic.detail import DetailView
from django.views.generic.list import ListView
from django.views.generic.detail import SingleObjectMixin

from ..forms import SubmissionForm, ContainerForm
from ..models import Submission


class OwnerMatchMixin(LoginRequiredMixin, UserPassesTestMixin):
    def test_func(self):
        obj = self.get_object()
        return obj.user == self.request.user


class SubmissionDetail(OwnerMatchMixin, DetailView):
    model = Submission
    # context_object_name = 'submission'
    DETAIL_FIELD_NAMES = (
        "ranked",
        "category",
        "url",
        "compute_time",
        "computing_and_hardware",
        "software",
        "method",
    )

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["container"] = self.object.container
        context["submission_details"] = [
            (field_name, getattr(self.object, field_name))
            for field_name in self.DETAIL_FIELD_NAMES
            if getattr(self.object, field_name)
        ]
        context["missing_fields"] = [
            field_name
            for field_name in self.DETAIL_FIELD_NAMES
            if not getattr(self.object, field_name)
        ]
        return context


class SubmissionList(LoginRequiredMixin, ListView):
    model = Submission

    def get_queryset(self):
        return super().get_queryset().filter(user=self.request.user)


class SubmissionDelete(OwnerMatchMixin, DeleteView):
    model = Submission
    success_url = reverse_lazy("submission-list")


@login_required
def edit_submission_view(request, pk=None):
    if request.method == "POST":
        # TODO: transaction?
        submission = Submission.objects.get(pk=pk) if pk else None
        container = submission.container if submission else None
        container_form = ContainerForm(request.POST, instance=container)
        submission_form = SubmissionForm(request.POST, instance=submission)
        if container_form.is_valid():
            container = container_form.save(commit=False)
            container.user = request.user
            container.save()
            if submission_form.is_valid():
                submission = submission_form.save(commit=False)
                submission.container = container
                submission.challenge = container.challenge
                submission.user = request.user
                submission.save()
                return redirect("submission-detail", pk=submission.pk)

    # context = {"container_form": container_form, "submission_form": submission_form}
    #  return render(request, "core/submission_form.html", context)

    elif request.method == "GET":
        if pk:
            submission = Submission.objects.get(pk=pk)
            container_form = ContainerForm(instance=submission.container)
            submission_form = SubmissionForm(instance=submission)
        else:
            container_form = ContainerForm()
            submission_form = SubmissionForm()
    else:
        return HttpResponseBadRequest()
    context = {"container_form": container_form, "submission_form": submission_form}
    return render(request, "core/submission_form.html", context)
