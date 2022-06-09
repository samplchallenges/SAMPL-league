# pylint: skip-file
import os
from pathlib import Path

from .base_settings import *  # lgtm [py/polluting-import]

DEBUG = False
# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent
MEDIA_ROOT = BASE_DIR / "media"
SECRET_KEY = "7&gz49!+qmxwltqo*!_g$n+i)%qcn9c^%2kbwzlnwiuofj29!%"

CRISPY_FAIL_SILENTLY = False

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql_psycopg2",
        "NAME": os.environ["RDS_DB_NAME"],
        "USER": os.environ["RDS_USERNAME"],
        "PASSWORD": os.environ["RDS_PASSWORD"],
        "HOST": os.environ["RDS_HOSTNAME"],
        "PORT": os.environ["RDS_PORT"],
    }
}
WORKER_CORES = int(os.environ.get("WORKER_CORES", 1))
WORKER_MEMORY = os.environ.get("WORKER_MEMORY", "4 GB")
WORKER_PROCESSES = int(os.environ.get("WORKER_PROCESSES", 1))
WORKER_QUEUE_PARTITION = os.environ.get("WORKER_QUEUE_PARTITION", "free")
WORKER_WALLTIME = os.environ.get("WORKER_WALLTIME", "12:00:00")

ECR_BASE_URL = os.environ["ECR_BASE_URL"]
ECR_SAMPLLEAGUE_URL = os.environ["ECR_SAMPLLEAGUE_URL"]
