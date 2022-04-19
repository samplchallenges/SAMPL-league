# pylint: skip-file
from .settings_prod import *  # lgtm [py/polluting-import]

BASE_DIR = Path("/data/homezvol0/osatom/")
MEDIA_ROOT = BASE_DIR / "sampl-app-extras/media"
LOGS_ROOT = BASE_DIR / "sampl-app-extras/logs"
SAMPL_ROOT = BASE_DIR / "SAMPL-league/"

# This settings file is on the dask workers, scheduler, and for the job that submits dask tasks
REMOTE_SCHEDULER = False


CONTAINER_ENGINE = "singularity"

# Dask SLURMCluster Settings
MINIMUM_WORKERS = 0
MAXIMUM_WORKERS = 1

# HPC3 Job Submitter Settings
JOB_SUBMITTER_LIFETIME = 60 * 60  # in seconds (1 hr)
CHECK_INTERVAL = 60  # in seconds (5 min)


LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "file": {
            "level": "DEBUG",
            "class": "logging.FileHandler",
            "filename": LOGS_ROOT / "app.log",
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
        "job_submitter": {
            "handlers": ["file", "rollbar"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "INFO"),
            "propagate": False,
        },
    },
}
