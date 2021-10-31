# pylint: disable=unused-argument, unused-variable

import pytest
from django.db import connection, reset_queries
from pytest_django.asserts import assertContains

from core.views import challenge as challenge_views


@pytest.mark.django_db
def test_challenge_detail(client, user, benzene_from_mol):
    client.force_login(user)
    challenge = benzene_from_mol.challenge
    response = client.get(f"/challenge/{challenge.pk}/")
    assertContains(response, "--molfile ChEBI_16716")


@pytest.mark.django_db
def test_challenge_list(client, user, benzene_from_mol):
    client.force_login(user)
    response = client.get("/challenge/")
    assertContains(response, "molfile_molw")


@pytest.mark.django_db
def test_query_count_elements_context(settings, benzene_from_mol):
    settings.DEBUG = True  # So that connection.queries will have data
    reset_queries()
    elements, input_types = challenge_views._input_elements(benzene_from_mol.challenge)
    query_count = len(connection.queries)
    expected_query_count = (
        1 + 1  # for list of input values  # for getting the file value
    )

    assert query_count == expected_query_count

    assert len(elements) == 1
    assert len(input_types) == 1
