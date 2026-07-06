"""
When we have files stored on S3, we want to seamlessly manage local copies for submission to containers.
"""

import logging
import os.path
import shutil

from django.conf import settings

BUFSIZE = 10 * 1024 * 1024  # Break big files into 10MB chunks
# BUFSIZE = 100

logger = logging.getLogger(__name__)


def _local_cache_path(field_file):
    local_path = os.path.join(settings.MEDIA_ROOT, field_file.name)
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    return local_path


def ensure_local_copy(field_file):
    """
    Should be a no-op if using default Django storage
    """
    local_path = _local_cache_path(field_file)
    if not os.path.exists(local_path):
        logger.debug("Creating local cache %s", local_path)
        with open(local_path, "wb") as local_fp:
            while True:
                data = field_file.read(BUFSIZE)
                if not data:
                    break
                local_fp.write(data)
                logger.debug("Writing chunk")

    return local_path


def preserve_local_copy(field_file, filepath):
    """
    If using local storage, no-op. If using remote storage, copy the file
    at filepath into the local cache path.
    """
    local_path = _local_cache_path(field_file)
    if not os.path.exists(local_path):
        logger.debug("Copying local file %s into %s", filepath, local_path)
        shutil.copyfile(filepath, local_path)


def delete_local_cache(field_file):
    local_path = _local_cache_path(field_file)
    if os.path.exists(local_path):
        # Only delete if using remote storage (e.g. S3!)
        try:
            _filepath = field_file.path
            logger.debug("Will not delete local copy as we don't have remote storage")
        except NotImplementedError:
            logger.debug("Deleting local copy as we have remote storage")
            os.remove(local_path)
