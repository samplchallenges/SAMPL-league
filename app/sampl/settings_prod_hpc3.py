# pylint: skip-file
import os

from .settings_prod import *  # lgtm [py/polluting-import]

BASE_DIR = Path(os.environ["BASE_DIR"])
MEDIA_ROOT = Path(os.environ["SAMPL_MEDIA_ROOT"])
LOGS_ROOT = Path(os.environ["SAMPL_LOGS_ROOT"])
DASK_WORKER_LOGS_ROOT = LOGS_ROOT / "dask-worker"
SAMPL_ROOT = Path(os.environ["SAMPL_ROOT"])  # /path/to/SAMPL-league

# This settings file is on the dask workers, scheduler, and for the job that submits dask tasks
REMOTE_SCHEDULER = False

# The container engine must be singularity on hpc3
CONTAINER_ENGINE = "singularity"

# HPC3 Job Submitter Settings
JOBQUEUE_CONFIG_FILE = SAMPL_ROOT / "app/referee/jobqueue.yaml"
JOB_SUBMITTER_WALLTIME = os.environ.get("JOB_SUBMITTER_WALLTIME", "24:00:00")
WALLTIME_HRS, WALLTIME_MIN, WALLTIME_SEC = tuple(JOB_SUBMITTER_WALLTIME.split(":"))
JOB_SUBMITTER_LIFETIME = (
    int(WALLTIME_HRS) * 60 * 60 + int(WALLTIME_MIN) * 60 + int(WALLTIME_SEC)
)
CHECK_INTERVAL = int(os.environ.get("CHECK_INTERVAL", 60))  # in seconds (1 min)

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
