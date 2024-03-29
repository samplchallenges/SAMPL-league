# pylint:disable=unused-argument
import os
import re
from datetime import datetime
from unittest.mock import Mock, patch

import dask.distributed as dd
import pytest
from django.core.management import call_command
from django.db import transaction
from django.forms.fields import CharField
from django.http.response import FileResponse
from django.urls import reverse
from django.utils import timezone

from core.forms import ContainerForm, SubmissionForm
from core.models import Challenge, ContainerType, Status, Submission
from core.tests import mocktime


@pytest.mark.django_db
def test_list_submissions(client, user, other_user, draft_submission):
    client.force_login(user)
    response = client.get("/submission/")
    assert response.status_code == 200
    assert draft_submission.name in response.content.decode()
    client.force_login(other_user)
    response = client.get("/submission/")
    assert response.status_code == 200
    assert draft_submission.name not in response.content.decode()


@patch("django.utils.timezone.now", mocktime.active)
@pytest.mark.django_db
def test_load_submission_form(client, user, draft_submission):
    client.force_login(user)
    response = client.get(
        f"/submission/add/?challenge_id={draft_submission.challenge_id}"
    )
    assert response.status_code == 200


@pytest.mark.django_db
def test_get_submission(client, user, other_user, draft_submission):
    client.force_login(user)
    response = client.get(f"/submission/{draft_submission.pk}/")
    assert response.status_code == 200
    assert response.context["submission"].draft_mode == True
    client.force_login(other_user)
    response = client.get(f"/submission/{draft_submission.pk}/")
    assert response.status_code == 403


@pytest.mark.django_db
def test_edit_submission(client, user, draft_submission):
    client.force_login(user)
    response = client.get(f"/submission/{draft_submission.pk}/edit/")
    submission_form = response.context["submission_form"]
    assert draft_submission.pk == submission_form.instance.pk


@pytest.mark.django_db
def test_clone_submission(client, user, draft_submission):
    client.force_login(user)
    response = client.get(f"/submission/{draft_submission.pk}/clone/")
    assert response.status_code == 200

    submission_form = response.context["submission_form"]

    assert submission_form.instance.id is None
    assert submission_form.instance.name == "Draft Submission"


@patch("django.utils.timezone.now", mocktime.active)
@pytest.mark.django_db
def test_create_submission(client, user, challenge):
    form_data = {
        "container-name": "test container",
        "container-container_type": ContainerType.DOCKER,
        "container-challenge": challenge.pk,
        "container-registry": "local",
        "container-label": "package2",
        "submission-name": "draft submission 2",
        "args-TOTAL_FORMS": 0,
        "args-INITIAL_FORMS": 0,
    }
    client.force_login(user)
    response = client.post("/submission/add/", form_data)
    assert response.status_code == 302
    response = client.get(response.url)
    submission = response.context["submission"]
    assert submission.draft_mode
    assert submission.name == form_data["submission-name"]


@patch("django.utils.timezone.now", mocktime.active)
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
    form_data.update({"args-TOTAL_FORMS": 0, "args-INITIAL_FORMS": 0})
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
    assert submission.draft_mode


@patch("django.utils.timezone.now", mocktime.inactive_after)
@pytest.mark.django_db
def test_load_expired_submission(client, user, draft_submission):
    client.force_login(user)
    response = client.get(f"/submission/{draft_submission.pk}/edit/")
    assert response.status_code == 200

    submission = response.context["submission"]
    assert not submission.challenge.is_active()

    # assert form submission, container, and args fields are visually disabled
    submission_form = response.context["submission_form"]
    container_form = response.context["submission_form"]
    arg_formset = response.context["arg_formset"]
    for field in submission_form.fields.keys():
        assert submission_form.fields[field].disabled
    for field in container_form.fields.keys():
        assert container_form.fields[field].disabled
    for form in arg_formset.forms:
        assert form.fields["key"].disabled
        assert form.fields["file_value"].disabled
        assert form.fields["DELETE"].disabled

    # assert submission note field is enabled
    submission_notes_form = response.context["submission_notes_form"]
    assert not submission_notes_form.fields["notes"].disabled


@patch("django.utils.timezone.now", mocktime.inactive_after)
@pytest.mark.django_db
def test_update_expired_submission(client, user, draft_submission):

    client.force_login(user)
    response = client.get(f"/submission/{draft_submission.pk}/edit/")
    assert response.status_code == 200

    submission = response.context["submission"]
    assert not submission.challenge.is_active()

    submission_form = response.context["submission_form"]
    container_form = response.context["submission_form"]
    submission_notes_form = response.context["submission_notes_form"]

    form_data = {f"{submission_notes_form.prefix}-notes": "hello from expired"}
    response = client.post(f"/submission/{draft_submission.pk}/edit/", form_data)
    response = client.get(f"/submission/{draft_submission.pk}/edit/")
    assert response.status_code == 200
    submission = response.context["submission"]
    assert submission.notes == "hello from expired"

    # update submission and container fields after expired - should not update
    submission_old = submission
    container_old = submission.container
    change_challenge = Challenge(
        name="ChangedChallenge",
        start_at=datetime(2020, 2, 1, hour=1, tzinfo=timezone.utc),
        end_at=datetime(2020, 9, 1, hour=1, tzinfo=timezone.utc),
        repo_url="http://github.com",
    )
    change_challenge.save()
    container_form_data = {
        "name": "UPDATED CONTAINER",
        "container_type": ContainerType.DOCKER,
        "registry": "docker.io",
        "challenge": change_challenge,
        "label": "osatom/adv-tutorial",
        "tag": "UPDATED",
    }
    submission_form_data = {
        "ranked": not submission_old.ranked,
        "url": "http://samplchallengesarecool.org",
        "method": "updating method from expired",
        "compute_time": "updating compute_time from expired",
        "software": "updating software from expired",
        "computing_and_hardware": "updating computing and hardware from expired",
    }
    final_container_form_data = {
        f"{container_form.prefix}-{key}": value
        for key, value in container_form_data.items()
    }
    final_submission_form_data = {
        f"{submission_form.prefix}-{key}": value
        for key, value in submission_form_data.items()
    }

    final_form_data = {**final_container_form_data, **final_submission_form_data}

    response = client.post(f"/submission/{draft_submission.pk}/edit/", final_form_data)
    assert response.status_code == 302

    response = client.get(f"/submission/{draft_submission.pk}/edit/")
    assert response.status_code == 200

    submission = response.context["submission"]
    container = submission.container

    for key in container_form_data.keys():
        assert getattr(container, key) == getattr(container_old, key)

    for key in submission_form_data.keys():
        assert getattr(submission, key) == getattr(submission_old, key)


@pytest.mark.django_db(transaction=True)
def test_run_submission(client, container_engine):

    # Because we have dask worker in a separate thread, we need to commit our transaction.
    # But the transaction test case will wipe out data from django's ContentTypes
    # So rerun our migrations to re-add our content types
    # Generally we don't run with processes False since it's a less thorough test
    # But for debugging with pdb, it's more convenient
    with patch("django.conf.settings.CONTAINER_ENGINE", container_engine):
        processes = True
        if processes:
            transaction.commit()
            call_command("migrate", "core", "zero", verbosity=0, interactive=False)
            call_command("migrate", "core", verbosity=0, interactive=False)
            call_command("sample_data")
            transaction.commit()
        else:
            call_command("sample_data")
        submission = Submission.objects.first()
        client.force_login(submission.user)

        future = None

        def save_future(l_future):
            nonlocal future
            future = l_future

        if processes:
            os.environ["CONTAINER_ENGINE"] = container_engine
            preload_file = "daskworkerinit_tst.py"
            cluster = dd.LocalCluster(n_workers=4, preload=(preload_file,))
        else:
            cluster = dd.LocalCluster(
                n_workers=1, processes=False, threads_per_worker=1
            )

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

        evaluation = public_run.evaluation_set.first()

        evaluation_url = reverse("evaluation-detail", kwargs={"pk": evaluation.pk})

        response = client.get(evaluation_url)
        assert response.context["evaluation"].pk == evaluation.pk

        log_url = reverse("evaluation-log-out", kwargs={"pk": evaluation.pk})
        response = client.get(log_url)
        assert response.context["log"] == evaluation.log_stdout

        errlog_url = reverse("evaluation-log-err", kwargs={"pk": evaluation.pk})
        response = client.get(errlog_url)
        assert response.context["log"] == evaluation.log_stderr


def test_download_container_arg_file(client, user, draft_submission, custom_file_arg):
    client.force_login(user)
    response = client.get(f"/download_arg/{custom_file_arg.id}/")
    assert response.status_code == 200


def test_download_input_file(client, user, benzene_from_mol):
    client.force_login(user)
    input_value = benzene_from_mol.inputvalue_set.first()
    response = client.get(f"/download_input/{input_value.pk}/")
    assert response.status_code == 200
    assert isinstance(response, FileResponse)
    assert re.match(
        'attachment; filename="ChEBI_16716.*.mdl"',
        response.headers["Content-Disposition"],
    )


@pytest.mark.django_db
def test_cancel_request(
    client, user, other_user, draft_submission, submission_run_factory
):
    submission_run = submission_run_factory(draft_submission)
    client.force_login(other_user)
    response = client.post(f"/submissionrun/{submission_run.pk}/cancel/")
    assert response.status_code == 404
    client.force_login(user)
    response = client.post(f"/submissionrun/{submission_run.pk}/cancel/")
    assert response.status_code == 302
    submission_run.refresh_from_db()
    assert submission_run.status == Status.CANCEL_PENDING


@pytest.mark.django_db
def test_cancel_request_remote(client, user, draft_submission, submission_run_factory):
    with patch("django.conf.settings.REMOTE_SCHEDULER", True):
        submission_run = submission_run_factory(draft_submission)
        submission_run.status = Status.PENDING_REMOTE
        submission_run.save(update_fields=["status"])
        client.force_login(user)
        response = client.post(f"/submissionrun/{submission_run.pk}/cancel/")
        assert response.status_code == 302
        submission_run.refresh_from_db()
        assert submission_run.status == Status.CANCELLED


@patch("django.conf.settings.REMOTE_SCHEDULER", True)
@pytest.mark.django_db(transaction=True)
def test_remote_scheduler(client):
    processes = True
    if processes:
        transaction.commit()
        call_command("migrate", "core", "zero", verbosity=0, interactive=False)
        call_command("migrate", "core", verbosity=0, interactive=False)
        call_command("sample_data")
        transaction.commit()
    else:
        call_command("sample_data")

    submission = Submission.objects.first()
    client.force_login(submission.user)
    response = client.post(f"/submission/{submission.pk}/submit/", {})
    assert response.url == reverse("submission-detail", kwargs={"pk": submission.pk})
    assert response.status_code == 302
    detail_url = response.url
    response = client.get(detail_url)
    public_run = response.context["public_run"]
    assert public_run.status == Status.PENDING_REMOTE
