# pylint: skip-file
import os
from pathlib import Path

from .base_settings import *  # lgtm [py/polluting-import]
from .debug_logging import LOGGING

DEBUG = False
# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent
MEDIA_ROOT = BASE_DIR / "media"
CONTAINER_FILES_ROOT = (
    Path(os.environ["CONTAINER_FILES_ROOT"])
    if os.environ.get("CONTAINER_FILES_ROOT", None)
    else None
)

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

ECR_BASE_URL = os.environ["ECR_BASE_URL"]
ECR_SAMPLLEAGUE_URL = os.environ["ECR_SAMPLLEAGUE_URL"]
