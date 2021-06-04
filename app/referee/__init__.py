import dask.distributed as dd
from django.conf import settings


def get_client():
    return dd.Client(settings.DASK_SCHEDULER_URL)
