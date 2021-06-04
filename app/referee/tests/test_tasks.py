import dask.distributed as dd
import pytest

from django.db import transaction
from django.core.management import call_command

from core.models import Submission
from referee import tasks


@pytest.mark.django_db(transaction=True)
def test_run_and_score_submission():
    # This test will fail if run after another transaction=True test
    # See workaround in tests/test_views.py:test_run_submission

    transaction.commit()
    call_command("migrate", "core", "zero", interactive=False)
    call_command("migrate", "core", interactive=False)
    call_command("sample_data")
    transaction.commit()

    submission = Submission.objects.first()
    cluster = dd.LocalCluster(n_workers=4, preload=("daskworkerinit_tst.py",))
    dask_client = dd.Client(cluster)

    print(submission.id, submission)
    future = tasks.run_and_score_submission(dask_client, submission)
    result = future.result()
    assert result
