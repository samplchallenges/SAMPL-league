# pylint: disable=unused-argument, unused-variable

import re
from pathlib import Path
from unittest.mock import patch

import pytest
from django.conf import settings

from core import models
from core.tests import mocktime

TEST_DATA_PATH = Path(__file__).parent.parent / "data"


def test_container(scoring_container):
    assert scoring_container.uri == "ghcr.io/megosato/score-coords:latest"

    scoring_container.tag = None

    assert scoring_container.uri == "ghcr.io/megosato/score-coords"


@pytest.mark.parametrize(
    ["mocknow", "compareval"],
    [
        (mocktime.inactive_before, False),
        (mocktime.active, True),
        (mocktime.inactive_after, False),
    ],
)
def test_challenge(challenge, mocknow, compareval):
    with patch("django.utils.timezone.now", mocknow):
        assert challenge.is_active() == compareval


def test_container_arg(draft_submission, custom_string_arg, custom_file_arg):
    container = draft_submission.container
    assert container.custom_args() == {"stringarg": "hello world"}
    filepath = str(
        Path(settings.MEDIA_ROOT)
        / "container_args/{}/{}/filearg/example.*.txt".format(
            container.user_id, container.id
        )
    )
    assert re.match(filepath, container.custom_file_args()["filearg"])


def test_cancel_requested(draft_submission, submission_run_factory):
    submission_run = submission_run_factory(draft_submission)
    submission_run.status = models.Status.CANCEL_PENDING
    submission_run.save()

    assert submission_run.check_cancel_requested()
