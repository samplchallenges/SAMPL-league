from django import forms
from django.core.exceptions import ValidationError


from .models import Submission, Container


class ContainerForm(forms.ModelForm):
    class Meta:
        model = Container
        fields = ["name", "challenge", "registry", "label", "tag", "digest"]
        widgets = {
            "challenge": forms.HiddenInput(),
        }


class SubmissionForm(forms.ModelForm):
    class Meta:
        model = Submission
        fields = [
            "draft_mode",
            "name",
            "challenge",
            "container",
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
            "challenge": forms.HiddenInput(),
            "container": forms.HiddenInput(),
        }
