from django import forms
from django.core.exceptions import ValidationError


from .models import Submission


class SubmissionForm(forms.ModelForm):
    class Meta:
        model = Submission
        fields = [
            "draft_mode",
            "name",
            "ranked",
            "url",
            "compute_time",
            "computing_and_hardware",
            "software",
            "category",
            "method",
        ]
