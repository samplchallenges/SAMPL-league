import pytest
from django.db import transaction

from referee import tasks


@pytest.mark.django_db(transaction=True)
def test_run_and_score_submission(dask_client, draft_submission, input_elements):
    # This test will fail if run after another transaction=True test
    # See workaround in tests/test_views.py:test_run_submission
    draft_submission.save()
    transaction.commit()
    print(draft_submission.id, draft_submission)
    future = tasks.run_and_score_submission(dask_client, draft_submission)
    result = future.result()
    assert result
