# pylint: skip-file
from io import StringIO

import pytest
from django.core.management import call_command

from core import models


@pytest.mark.django_db
def test_sample_data():
    out = StringIO()
    call_command("sample_data", stdout=out, stderr=StringIO())

    challenge = models.Challenge.objects.get(name="Milo")
    assert challenge.container_set.count() == 2
    assert challenge.submission_set.count() == 1
    submission = models.Submission.objects.get()
    assert submission.submissionrun_set.count() == 1
    assert models.InputElement.objects.all().count() == 6
    assert models.ValueType.objects.all().count() == 2
    assert models.InputValue.objects.all().count() == 6
    assert models.FloatValue.objects.all().count() == 6
