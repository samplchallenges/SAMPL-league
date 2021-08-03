import os

from .settings_dev import *  # lgtm [py/polluting-import]

DEFAULT_FILE_STORAGE = "storages.backends.s3boto3.S3Boto3Storage"
AWS_ACCESS_KEY_ID = os.environ["AWS_ACCESS_KEY_ID_S3"]
AWS_SECRET_ACCESS_KEY = os.environ["AWS_SECRET_ACCESS_KEY_S3"]
AWS_STORAGE_BUCKET_NAME = os.environ["AWS_STORAGE_BUCKET_NAME"]
AWS_S3_REGION_NAME = "us-east-2"


LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
        },
    },
    "root": {
        "handlers": ["console"],
        "level": "INFO",
    },
    "loggers": {
        "django": {
            "handlers": ["console"],
            "level": os.getenv("DJANGO_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
        "core": {
            "handlers": ["console"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
        "referee": {
            "handlers": ["console"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
        "ever_given": {
            "handlers": ["console"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
    },
}
