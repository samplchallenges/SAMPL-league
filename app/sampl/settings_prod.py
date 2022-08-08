# pylint: skip-file
import os
from functools import partial
from pathlib import Path

from . import aws_login
from .base_settings import *  # lgtm [py/polluting-import]

BASE_DIR = Path("/opt/app/sampl")
MEDIA_ROOT = BASE_DIR / "media"

# Set dask URL for aws
# DASK_SCHEDULER_URL = "172.31.39.184:8786"

DASK_SCHEDULER_HOST = os.environ.get("DASK_SCHEDULER_HOST", "127.0.0.1")
DASK_SCHEDULER_URL = f"{DASK_SCHEDULER_HOST}:8786"


# We only need rollbar in production
MIDDLEWARE.append("rollbar.contrib.django.middleware.RollbarNotifierMiddleware")

DEBUG = False
SECRET_KEY = os.environ["SAMPL_SECRET_KEY"]
# Rollbar API key
POST_SERVER_ITEM_ACCESS_TOKEN = os.environ["POST_SERVER_ITEM_ACCESS_TOKEN"]

MAIN_WEBSITE_URL = os.environ.get("MAIN_WEBSITE_URL", "app.samplchallenges.org")

ALLOWED_HOSTS = [
    "127.0.0.1",
    "localhost",
    MAIN_WEBSITE_URL,
    os.environ.get("EC2_PRIVATE_IP", "172.31.43.228"),
    os.environ.get("SAMPL_DNS_1", ""),
    os.environ.get("SAMPL_DNS_2", ""),
]
CSRF_TRUSTED_ORIGINS = [
    f"https://{MAIN_WEBSITE_URL}",
]

DEFAULT_FILE_STORAGE = "storages.backends.s3boto3.S3Boto3Storage"
AWS_ACCESS_KEY_ID = os.environ["AWS_ACCESS_KEY_ID_S3"]
AWS_SECRET_ACCESS_KEY = os.environ["AWS_SECRET_ACCESS_KEY_S3"]
AWS_STORAGE_BUCKET_NAME = "sampl-league-storage"
AWS_S3_REGION_NAME = "us-east-2"

LOGIN_TO_AWS = True
ECR_BASE_URL = os.environ["ECR_BASE_URL"]
ECR_SAMPLLEAGUE_URL = os.environ["ECR_SAMPLLEAGUE_URL"]

AWS_LOGIN_FUNCTION = partial(aws_login.run_aws_login, ECR_BASE_URL)

AWS_LOGOUT_FUNCTION = aws_login.run_aws_logout

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

# rollbar settings
ROLLBAR = {
    "access_token": POST_SERVER_ITEM_ACCESS_TOKEN,
    "environment": "production",
    "branch": "main",
    "root": str(BASE_DIR),
}

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "file": {
            "level": "DEBUG",
            "class": "logging.FileHandler",
            "filename": "/var/log/app.log",
        },
        "rollbar": {
            "level": "DEBUG",
            "access_token": POST_SERVER_ITEM_ACCESS_TOKEN,
            "environment": "production",
            "class": "rollbar.logger.RollbarHandler",
        },
    },
    "loggers": {
        "django": {
            "handlers": ["file", "rollbar"],
            "level": os.getenv("DJANGO_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
        "core": {
            "handlers": ["file", "rollbar"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
        "referee": {
            "handlers": ["file", "rollbar"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
        "ever_given": {
            "handlers": ["file", "rollbar"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
    },
}


# which container engine to use, must be "docker" or "singularity"
CONTAINER_ENGINE = os.environ.get("CONTAINER_ENGINE", "docker")
