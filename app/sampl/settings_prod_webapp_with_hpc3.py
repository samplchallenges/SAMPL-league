# pylint: skip-file
from .settings_prod import *  # lgtm [py/polluting-import]

BASE_DIR = Path("/opt/app/sampl")
MEDIA_ROOT = BASE_DIR / "media"

# Settings file for webapp (running on AWS) when the
# dask workers and scheduler will be on HPC3
REMOTE_SCHEDULER = True

DEBUG = True

WORKER_CORES = int(os.environ.get("WORKER_CORES", 1))
WORKER_MEMORY = os.environ.get("WORKER_MEMORY", "4 GB")
WORKER_PROCESSES = int(os.environ.get("WORKER_PROCESSES", 1))
WORKER_QUEUE_PARTITION = os.environ.get("WORKER_QUEUE_PARTITION", "free")
WORKER_WALLTIME = os.environ.get("WORKER_WALLTIME", "12:00:00")
