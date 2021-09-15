from crispy_forms import layout
from crispy_forms.helper import FormHelper
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.forms import inlineformset_factory

from .models import Container, Submission, SubmissionArg


class RegisterForm(UserCreationForm):
    class Meta(UserCreationForm.Meta):
        fields = UserCreationForm.Meta.fields + ("email",)


class ContainerForm(forms.ModelForm):
    prefix = "container"

    class Meta:
        model = Container
        fields = ["name", "challenge", "registry", "label", "tag"]


class SubmissionForm(forms.ModelForm):
    prefix = "submission"

    class Meta:
        model = Submission
        fields = [
            "draft_mode",
            "name",
            "ranked",
            "category",
            "url",
            "compute_time",
            "computing_and_hardware",
            "software",
            "method",
        ]
        widgets = {
            "compute_time": forms.Textarea(attrs={"cols": 30, "rows": 4}),
            "software": forms.Textarea(attrs={"cols": 30, "rows": 4}),
        }


class SubmissionArgStringForm(forms.ModelForm):
    prefix = "arg_string"

    class Meta:
        model = SubmissionArg
        fields = ["key", "string_value"]


class SubmissionArgFileForm(forms.ModelForm):
    prefix = "arg_file"

    class Meta:
        model = SubmissionArg
        fields = ["key", "file_value"]


def submission_arg_string_formset():
    return inlineformset_factory(
        Submission,
        SubmissionArg,
        form=SubmissionArgStringForm,
        fields=("key", "string_value"),
        extra=1,
    )


def submission_arg_file_formset():
    return inlineformset_factory(
        Submission,
        SubmissionArg,
        form=SubmissionArgFileForm,
        fields=("key", "file_value"),
        extra=1,
    )


def submission_arg_formset():
    return inlineformset_factory(
        Submission, SubmissionArg, fields=("key", "file_value"), extra=1
    )


class ArgFormHelper(FormHelper):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.disable_csrf = True
        self.form_tag = False
        self.field_template = "core/arg_field.html"
        self.template = "core/arg_formset.html"
