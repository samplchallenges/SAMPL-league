from crispy_forms.helper import FormHelper
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.forms import inlineformset_factory
from django.utils import timezone
from .models import Container, ContainerArg, Submission
from . import configurator

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

class SubmissionNotesForm(forms.Form):
    prefix = "submission_note"
    notes = forms.CharField(label='Notes', widget=forms.Textarea(attrs={"cols": 30, "rows": 4, "placeholder": configurator.NOTES_DETAILS}), required=False)

def container_arg_formset():
    return inlineformset_factory(
        Container, ContainerArg, fields=("key", "file_value"), extra=1
    )


class ArgFormHelper(FormHelper):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.disable_csrf = True
        self.form_tag = False
        self.field_template = "core/arg_field.html"
        self.template = "core/arg_formset.html"
