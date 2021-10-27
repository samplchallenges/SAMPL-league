# pylint: disable=unused-argument, unused-variable
import re

import pytest
from pytest_django.asserts import assertContains


@pytest.mark.django_db
def test_challenge_detail(client, user, benzene_from_mol):
    client.force_login(user)
    challenge = benzene_from_mol.challenge
    response = client.get(f"/challenge/{challenge.pk}/")

    response = "".join(i for i in response.content.decode('utf-8'))
    regex = re.compile(r'--molfile ChEBI_16716.*.mdl')
    
    assert regex.search(response)

@pytest.mark.django_db
def test_challenge_list(client, user, benzene_from_mol):
    client.force_login(user)
    response = client.get("/challenge/")
    assertContains(response, "molfile_molw")
