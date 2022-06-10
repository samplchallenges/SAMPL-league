from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.http import HttpResponseBadRequest
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse_lazy
from django.views.generic.detail import DetailView
from django.views.generic.edit import DeleteView
from django.views.generic.list import ListView

import referee
import referee.tasks

from .. import forms
from ..models import Challenge, Submission, SubmissionRun

# pylint: disable=too-many-ancestors


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


class SubmissionDetail(OwnerMatchMixin, DetailView):
    model = Submission

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["container"] = self.object.container
        context["custom_args"] = self.object.container.args.all()
        context["submission_details"] = self.object.details
        context["missing_fields"] = self.object.missing_fields
        context["public_run"] = self.object.last_public_run()
        if context["public_run"]:
            context["public_completion"] = context["public_run"].completion()
        context["private_run"] = self.object.last_private_run()
        if context["private_run"]:
            context["private_completion"] = context["private_run"].completion()
        return context


class SubmissionList(LoginRequiredMixin, ListView):
    model = Submission

    def get_queryset(self):
        return super().get_queryset().filter(user=self.request.user)


class SubmissionDelete(OwnerMatchMixin, DeleteView):
    model = Submission
    success_url = reverse_lazy("submission-list")


def ignore_future(future):
    """patch this in testing"""
    return future.key


@login_required
def cancel_submissionrun_view(request, pk):
    if request.method != "POST":
        return HttpResponseBadRequest()
    submission_run = get_object_or_404(
        SubmissionRun, pk=pk, submission__user=request.user
    )
    if not submission_run.is_finished():
        submission_run.mark_for_cancel()

    return redirect("submission-detail", pk=submission_run.submission.pk)


@login_required
def submit_submission_view(request, pk):
    if request.method != "POST":
        return HttpResponseBadRequest()
    submission = get_object_or_404(Submission, pk=pk, user=request.user)

    if settings.REMOTE_SCHEDULER:
        referee.tasks.enqueue_submission(submission)
    else:
        # verifies that user matches
        dask_client = referee.get_client()
        future = referee.tasks.run_and_score_submission(dask_client, submission)
        ignore_future(future)
    return redirect("submission-detail", pk=submission.pk)


def _create_submission(container, submission_form, submission_notes_form):
    submission = submission_form.save(commit=False)
    submission.container = container
    submission.challenge = container.challenge
    submission.user = container.user
    if submission_notes_form.is_valid():
        submission.notes = submission_notes_form.cleaned_data["notes"]
    submission.save()
    return submission


def _create_container_args(container, arg_formset):
    arg_instances = arg_formset.save(commit=False)
    for instance in arg_instances:
        instance.container = container
        instance.save()
    for instance in arg_formset.deleted_objects:
        instance.delete()


def _disable_fields_for_inactive_challenge(
    submission_form, container_form, arg_formset
):
    for field in submission_form.fields.keys():
        if field != "notes":
            submission_form.fields[field].disabled = True

    for field in container_form.fields.keys():
        container_form.fields[field].disabled = True

    for form in arg_formset.forms:
        form.fields["key"].required = False
        form.fields["key"].disabled = True
        form.fields["file_value"].disabled = True
        form.fields["DELETE"].disabled = True


@login_required
def edit_submission_view(request, pk=None, clone=False):
    form_action = ""
    show_container = True
    show_args = False
    if request.method == "POST":
        # TODO: transaction?
        submission = Submission.objects.get(pk=pk, user=request.user) if pk else None
        if submission:
            challenge = submission.challenge
        container = submission.container if submission else None
        container_form = forms.ContainerForm(request.POST, instance=container)
        submission_form = forms.SubmissionForm(request.POST, instance=submission)
        submission_notes_form = forms.SubmissionNotesForm(request.POST)
        ArgFormSet = forms.container_arg_formset()
        arg_formset = None

        if pk:
            if not submission.challenge.is_active():
                if submission_notes_form.is_valid():
                    submission.notes = submission_notes_form.cleaned_data["notes"]
                    submission.save(update_fields=["notes"])
                return redirect("submission-detail", pk=submission.pk)
        if container_form.is_valid():
            challenge = container_form.cleaned_data["challenge"]
            if challenge.is_active():
                container = container_form.save(commit=False)
                container.user = request.user
                container.save()
                if submission_form.is_valid():
                    submission = _create_submission(
                        container, submission_form, submission_notes_form
                    )
                    arg_formset = ArgFormSet(
                        request.POST, request.FILES, instance=container
                    )
                    if arg_formset.is_valid():
                        _create_container_args(container, arg_formset)
                        return redirect("submission-detail", pk=submission.pk)
                    else:
                        show_args = True

        if arg_formset is None:
            arg_formset = ArgFormSet(request.POST, request.FILES, instance=container)
    elif request.method == "GET":
        if pk:
            show_container = False
            submission = Submission.objects.get(pk=pk)
            container = submission.container
            challenge = container.challenge
            if clone:
                submission.pk = None
                container.pk = None
                form_action = reverse_lazy("submission-add")
            container_form = forms.ContainerForm(instance=container)

            submission_form = forms.SubmissionForm(instance=submission)
            submission_notes_form = forms.SubmissionNotesForm(
                initial={"notes": submission.notes}
            )
            arg_formset = forms.container_arg_formset()(instance=container)

        else:
            initial_values = {}
            if "challenge_id" in request.GET:
                initial_values["challenge"] = request.GET["challenge_id"]
                challenge = Challenge.objects.get(pk=initial_values["challenge"])
                if not challenge.is_active():
                    return redirect("challenge-detail", pk=challenge.pk)
            else:
                return HttpResponseBadRequest()

            container_form = forms.ContainerForm(initial=initial_values)
            submission_form = forms.SubmissionForm()
            submission_notes_form = forms.SubmissionNotesForm()
            arg_formset = forms.container_arg_formset()()
            submission = None

        if pk and not submission.challenge.is_active():
            _disable_fields_for_inactive_challenge(
                submission_form, container_form, arg_formset
            )

    else:
        return HttpResponseBadRequest()
    context = {
        "container_form": container_form,
        "show_container": show_container,
        "show_args": show_args,
        "submission_form": submission_form,
        "submission_notes_form": submission_notes_form,
        "arg_formset": arg_formset,
        "arg_helper": forms.ArgFormHelper(),
        "form_action": form_action,
        "challenge": challenge,
        "submission": submission,
    }
    return render(request, "core/submission_form.html", context)
