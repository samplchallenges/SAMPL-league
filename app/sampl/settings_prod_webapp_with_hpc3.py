# pylint: skip-file
from .settings_prod import *  # lgtm [py/polluting-import]

BASE_DIR = Path("/opt/app/sampl")
MEDIA_ROOT = BASE_DIR / "media"

# Settings file for webapp (running on AWS) when the
# dask workers and scheduler will be on HPC3
REMOTE_SCHEDULER = True

DEBUG = True
