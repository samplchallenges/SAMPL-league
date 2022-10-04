# pylint: disable=unused-argument, unused-variable

import re
from pathlib import Path
from unittest.mock import patch

import pytest
from django.conf import settings

from core import models
from core.models.admin_managed import NotFullyLoadedException
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


def test_fully_loaded_elements(benzene_from_mol):
    benzene_from_mol.fully_loaded()

    aks = list(benzene_from_mol.answerkey_set.all())
    aks[0].delete()

    with pytest.raises(NotFullyLoadedException):
        benzene_from_mol.fully_loaded()


def test_fully_loaded_parent_element(smiles_docking_config_and_func):
    config, add_element = smiles_docking_config_and_func

    el1 = add_element("Benzene", "c1ccccc1", 78.01)
    el2 = add_element("Methane", "C", 16.04)

    parent_element = el1.parent

    parent_element.fully_loaded()
    with pytest.raises(NotFullyLoadedException):
        # Not fully loaded until we have private elements
        config.challenge.fully_loaded()

    el2.is_public = False
    el2.save()
    config.challenge.fully_loaded()

    aks = list(el1.answerkey_set.all())
    aks[0].delete()

    with pytest.raises(NotFullyLoadedException):
        parent_element.fully_loaded()

    with pytest.raises(NotFullyLoadedException):
        config.challenge.fully_loaded()


def test_fully_loaded_challenge_scoremaker(smiles_molw_config, input_elements):
    challenge = smiles_molw_config.challenge
    challenge.fully_loaded()

    challenge.scoremaker.delete()
    challenge.refresh_from_db()
    with pytest.raises(NotFullyLoadedException):
        challenge.fully_loaded()


def test_fully_loaded_challenge_scoretypes(smiles_molw_config, input_elements):
    challenge = smiles_molw_config.challenge
    challenge.fully_loaded()

    challenge.scoretype_set.first().delete()
    challenge = models.Challenge.objects.get(
        id=challenge.id
    )  # Clear cached_property value

    with pytest.raises(NotFullyLoadedException):
        challenge.fully_loaded()
