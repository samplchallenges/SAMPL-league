from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.http import HttpResponseBadRequest
from django.shortcuts import redirect, render
from django.urls import reverse_lazy
from django.views.generic.detail import DetailView
from django.views.generic.edit import DeleteView
from django.views.generic.list import ListView

import referee.tasks as referee

from ..forms import ContainerForm, SubmissionForm
from ..models import Submission


class OwnerMatchMixin(LoginRequiredMixin, UserPassesTestMixin):
    """
    This mixin on a view will make sure the user accessing the object in the view
    is the owner of that object. test_func is the Django builtin that will be used
    in the standard Django mixin 'UserPassesTestMixin'
    """

    def test_func(self):
        if self.request.user.is_superuser:
            return True
        obj = self.get_object()
        return obj.user == self.request.user


# pylint: disable=too-many-ancestors


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
        try:
            context["public_run"] = self.object.submissionrun_set.filter(
                is_public=True
            ).latest("updated_at")
            context["public_status"] = referee.get_status(context["public_run"].key)
            context["private_run"] = self.object.submissionrun_set.filter(
                is_public=False
            ).latest("updated_at")
        except:
            pass
        return context


class SubmissionList(LoginRequiredMixin, ListView):
    model = Submission

    def get_queryset(self):
        return super().get_queryset().filter(user=self.request.user)


class SubmissionDelete(OwnerMatchMixin, DeleteView):
    model = Submission
    success_url = reverse_lazy("submission-list")


@login_required
def submit_submission_view(request, pk):
    if request.method != "POST":
        return HttpResponseBadRequest()
    submission = Submission.objects.get(pk=pk, user=request.user)
    # verifies that user matches
    submission_run_1 = referee.run_submission(submission.pk, is_public=True)
    submission_run_2 = referee.run_submission(submission.pk, is_public=False)
    referee.score_submission(submission.pk, submission_run_1.pk, submission_run_2.pk)
    return redirect("submission-detail", pk=submission.pk)


@login_required
def clone_submission_view(request, pk):
    if request.method == "GET":
        submission = Submission.objects.get(pk=pk, user=request.user).clone()
        submission.save()
        return redirect("submission-update", pk=submission.pk)

    return HttpResponseBadRequest()


@login_required
def edit_submission_view(request, pk=None, clone=False):
    form_action = ""
    if request.method == "POST":
        # TODO: transaction?
        submission = Submission.objects.get(pk=pk, user=request.user) if pk else None
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
    elif request.method == "GET":
        if pk:
            submission = Submission.objects.get(pk=pk)
            if clone:
                submission.pk = None
                submission.container.pk = None
                form_action = reverse_lazy("submission-add")
            container_form = ContainerForm(instance=submission.container)
            submission_form = SubmissionForm(instance=submission)
        else:
            container_form = ContainerForm()
            submission_form = SubmissionForm()
    else:
        return HttpResponseBadRequest()
    context = {
        "container_form": container_form,
        "submission_form": submission_form,
        "form_action": form_action,
    }
    return render(request, "core/submission_form.html", context)
