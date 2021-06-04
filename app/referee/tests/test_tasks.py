import pytest
from django.db import transaction
from django.core.management import call_command

from referee import tasks


@pytest.mark.django_db(transaction=True)
def test_run_and_score_submission(dask_client, draft_submission, input_elements):
    # This test will fail if run after another transaction=True test
    # See workaround in tests/test_views.py:test_run_submission

    transaction.commit()
    call_command("migrate", "core", "zero", interactive=False)
    call_command("migrate", "core", interactive=False)
    call_command("sample_data")
    transaction.commit()

    print(draft_submission.id, draft_submission)
    future = tasks.run_and_score_submission(dask_client, draft_submission)
    result = future.result()
    assert result
