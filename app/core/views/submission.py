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


class SubmissionList(LoginRequiredMixin, ListView):
    model = Submission

    def get_queryset(self):
        return super().get_queryset().filter(user=self.request.user)


class SubmissionDelete(OwnerMatchMixin, DeleteView):
    model = Submission
    success_url = reverse_lazy("submission-list")


@login_required
def edit_submission_view(request):
    if request.method == "POST":
        # TODO: transaction?
        container_form = ContainerForm(request.POST)
        submission_form = SubmissionForm(request.POST)
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

        # raise Exception(
        #    "A form was invalid: container  "
        #    + str(container_form.errors)
        #    + "\n: submission "
        #    + str(submission_form.errors)
        # )
        context = {"container_form": container_form, "submission_form": submission_form}
        return render(request, "core/submission_form.html", context)

    elif request.method == "GET":
        container_form = ContainerForm()
        submission_form = SubmissionForm()
    else:
        return HttpResponseBadRequest()
    context = {"container_form": container_form, "submission_form": submission_form}
    return render(request, "core/submission_form.html", context)
