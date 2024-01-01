# pylint: disable=redefined-outer-name

import pytest
from django.contrib.auth import get_user_model
from django.test import Client

# from django.urls import reverse


@pytest.fixture
def staff_user(db):  # pylint: disable=unused-argument
    User = get_user_model()
    user = User.objects.create_superuser(username="staffer", password="hello")
    return user


@pytest.fixture
def c(staff_user):
    c = Client()
    c.force_login(staff_user)
    return c


def test_challenge_admin(molfile_molw_config, c):
    response = c.get(
        f"/admin/core/challenge/{molfile_molw_config.challenge.id}/change/"
    )
    assert response.status_code == 200

    assert b"No batches yet" in response.content

    molfile_molw_config.challenge.max_batch_size = 3
    molfile_molw_config.challenge.save()

    response = c.get(
        f"/admin/core/challenge/{molfile_molw_config.challenge.id}/change/"
    )
    assert response.status_code == 200

    assert b"batches, created" in response.content


def test_input_value(benzene_from_mol, c):
    response = c.get(f"/admin/core/inputvalue/{benzene_from_mol.id}/change/")
    assert response.status_code == 200

    html = response.content.decode("utf8")

    assert "Object link" in html
    assert "Value type challenge" in html
