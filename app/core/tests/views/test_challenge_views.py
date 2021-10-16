# pylint: disable=unused-argument, unused-variable
import re
import pytest
from pytest_django.asserts import assertContains


@pytest.mark.django_db
def test_challenge_detail(client, user, benzene_from_mol):
    client.force_login(user)
    challenge = benzene_from_mol.challenge
    response = client.get(f"/challenge/{challenge.pk}/")

    found = False
    for elem in response.context[0].dicts[3]['elements'][0]:
        print(elem, type(elem))

        if re.match("--molfile ChEBI_16716.*.mdl", elem.strip()):
            found = True
            break
    assert found


@pytest.mark.django_db
def test_challenge_list(client, user, benzene_from_mol):
    client.force_login(user)
    response = client.get("/challenge/")
    assertContains(response, "molfile_molw")
