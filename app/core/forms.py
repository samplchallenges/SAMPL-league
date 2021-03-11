from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.core.exceptions import ValidationError


from .models import Submission, Container


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
