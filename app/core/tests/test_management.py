# pylint: disable=unused-argument, unused-variable
from io import StringIO
from pathlib import Path

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


@pytest.mark.django_db
def test_load_yaml():
    config_data_dir = Path(__file__).parent.parent.parent / "config_data"
    config_yaml = config_data_dir / "docking" / "docking.challenge.yml"
    out = StringIO()
    call_command("load_yaml", config_yaml, stdout=out, stderr=StringIO())
    challenge = models.Challenge.objects.get(name="Example Docking Challenge")
    challenge.fully_loaded()

    assert challenge.inputelement_set.filter(is_parent=True).count() == 1
    for is_public in (True, False):
        assert (
            challenge.inputelement_set.filter(
                is_parent=False, is_public=is_public
            ).count()
            == 3
        )

    with pytest.raises(Exception):
        call_command("load_yaml", config_yaml, stdout=out, stderr=StringIO())
    call_command("load_yaml", config_yaml, "--delete", stdout=out, stderr=StringIO())
