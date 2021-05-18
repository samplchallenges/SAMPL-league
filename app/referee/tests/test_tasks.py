import pytest
import dask.distributed as dd
from django.db import transaction

from referee import tasks


@pytest.fixture(scope="session")
def dask_client():
    # We need at least 4 workers
    with dd.LocalCluster(n_workers=4, preload=("daskworkerinit_tst.py",)) as cluster:
        yield dd.Client(cluster)


@pytest.mark.django_db(transaction=True)
def test_run_and_score_submission(dask_client, draft_submission):
    draft_submission.save()
    transaction.commit()
    print(draft_submission.id, draft_submission)
    future = tasks.run_and_score_submission(dask_client, draft_submission)
    result = future.result
    assert False
