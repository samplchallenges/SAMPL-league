import dask.distributed as dd
from django.conf import settings


def get_client():
    """For interactive troubleshooting:
    cluster = dd.LocalCluster(n_workers=1, processes=False, threads_per_worker=1)

    return dd.Client(cluster)
    """
    return dd.Client(settings.DASK_SCHEDULER_URL)
