from .settings_prod import *  # lgtm [py/polluting-import]

BASE_DIR = Path("/opt/app/sampl")
MEDIA_ROOT = BASE_DIR / "media"

# This settings file is on the dask workers, scheduler, and for the job that submits dask tasks
REMOTE_SCHEDULER = False
