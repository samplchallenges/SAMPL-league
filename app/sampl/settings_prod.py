import os

from .base_settings import *  # lgtm [py/polluting-import]

# Set dask URL for aws
DASK_SCHEDULER_URL = "172.31.39.184:8786"

# We only need rollbar in production
MIDDLEWARE.append("rollbar.contrib.django.middleware.RollbarNotifierMiddleware")

DEBUG = False
SECRET_KEY = os.environ["SAMPL_SECRET_KEY"]
# Rollbar API key
POST_SERVER_ITEM_ACCESS_TOKEN = os.environ["POST_SERVER_ITEM_ACCESS_TOKEN"]
ALLOWED_HOSTS = [
    "app.samplchallenges.org",
    "sampl.us-east-2.elasticbeanstalk.com",
]

DEFAULT_FILE_STORAGE = "storages.backends.s3boto3.S3Boto3Storage"
AWS_ACCESS_KEY_ID = os.environ["AWS_ACCESS_KEY_ID_S3"]
AWS_SECRET_ACCESS_KEY = os.environ["AWS_SECRET_ACCESS_KEY_S3"]
AWS_STORAGE_BUCKET_NAME = "sampl-league-storage"
AWS_S3_REGION_NAME = "us-east-2"

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
    "root": BASE_DIR,
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
