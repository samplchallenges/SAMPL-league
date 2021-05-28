from unittest.mock import Mock, patch

import pytest

import dask.distributed as dd

from django.core.management import call_command
from django.db import transaction
from django.forms.fields import CharField
from django.urls import reverse

from core.forms import ContainerForm, SubmissionForm
from core.views.submission import edit_submission_view
from core.models import Submission


@pytest.mark.django_db
def test_load_submission_form(rf, user, draft_submission):
    request = rf.get("/core/submission/add/")
    request.user = user
    response = edit_submission_view(request)
    assert response.status_code == 200


@pytest.mark.django_db
def test_get_submission(rf, user, draft_submission):
    request = rf.get(f"/core/submission/{draft_submission.pk}/")
    request.user = user
    response = edit_submission_view(request)
    assert response.status_code == 200


@pytest.mark.django_db
def test_create_submission(client, user, challenge):
    form_data = {
        "container-name": "test container",
        "container-challenge": challenge.pk,
        "container-registry": "local",
        "container-label": "package2",
        "submission-name": "draft submission 2",
    }
    client.force_login(user)
    response = client.post("/submission/add/", form_data)

    assert response.status_code == 302
    response = client.get(response.url)
    submission = response.context["submission"]
    assert submission.draft_mode
    assert submission.name == form_data["submission-name"]


@pytest.mark.django_db
def test_update_submission(client, user, draft_submission):
    client.force_login(user)
    submission_form = SubmissionForm(instance=draft_submission)
    container_form = ContainerForm(instance=draft_submission.container)

    def skey(key):
        return f"{submission_form.prefix}-{key}"

    def ckey(key):
        return f"{container_form.prefix}-{key}"

    form_data = {
        skey(key): "sample data"
        for key, field in submission_form.fields.items()
        if isinstance(field, CharField)
    }
    form_data[skey("url")] = "http://test.org"
    form_data[skey("draft_mode")] = False
    form_data[skey("category")] = "Mixed"
    form_data[skey("ranked")] = False
    form_data.update(
        {skey(key): value for key, value in submission_form.initial.items() if value}
    )
    response = client.post(f"/submission/{draft_submission.pk}/edit/", form_data)
    assert response.status_code == 200
    assert response.context["submission_form"].is_valid()
    assert not response.context["container_form"].is_valid()
    form_data.update(
        {ckey(key): value for key, value in container_form.initial.items() if value}
    )
    response = client.post(f"/submission/{draft_submission.pk}/edit/", form_data)
    assert response.status_code == 302
    assert response.url == reverse(
        "submission-detail", kwargs={"pk": draft_submission.pk}
    )
    response = client.get(response.url)
    submission = response.context["submission"]
    assert not submission.draft_mode


@pytest.mark.django_db(transaction=True)
def test_run_submission(client):

    # Because we have dask worker in a separate thread, we need to commit out transaction.
    # But the transaction test case will wipe out data from django's ContentTypes
    # So rerun our migrations to re-add our content types

    transaction.commit()
    call_command("migrate", "core", "zero", interactive=False)
    call_command("migrate", "core", interactive=False)
    call_command("sample_data")
    transaction.commit()

    submission = Submission.objects.first()
    client.force_login(submission.user)

    future = None

    def save_future(l_future):
        nonlocal future
        future = l_future

    cluster = dd.LocalCluster(n_workers=4, preload=("daskworkerinit_tst.py",))
    dask_client = dd.Client(cluster)

    mock_get_client = Mock(return_value=dask_client)
    with patch("referee.get_client", mock_get_client):
        with patch("core.views.submission.ignore_future", save_future):
            response = client.post(f"/submission/{submission.pk}/submit/", {})

        assert response.status_code == 302
        assert response.url == reverse(
            "submission-detail", kwargs={"pk": submission.pk}
        )
        detail_url = response.url
        response = client.get(detail_url)

        # public_run = response.context["public_run"]
        # if public_run is not None:
        #   assert public_run.status in ("PENDING", "SUCCESS")
        result = future.result()
        response = client.get(detail_url)
        public_run = response.context["public_run"]
        assert public_run.status == "SUCCESS"
        assert result

    cluster.close()
