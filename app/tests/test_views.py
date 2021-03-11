import pytest

from django.forms.fields import CharField
from django.urls import reverse

from core.forms import ContainerForm, SubmissionForm
from core.views.submission import edit_submission_view


@pytest.mark.django_db
def test_load_submission_form(rf, user, draft_submission):
    request = rf.get(f"/core/submission/add/")
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
