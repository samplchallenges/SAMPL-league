"""
When we have files stored on S3, we want to seamlessly manage local copies for submission to containers.
"""
import logging
import os.path

import boto3


BUFSIZE = 10*1024*1024  # Break big files into 10MB chunks


LOGGER = logging.getLogger(__name__)


def ensure_local_copy(field_file):
    """
    Should be a no-op if using default Django storage
    """
    local_path = field_file.path
    if not os.path.exists(local_path):
        LOGGER.debug("Creating local cache %s", local_path)
        with open(local_path, 'wb') as local_fp:
            while True:
                data = field_file.read(BUFSIZE)
                if not data:
                    break
                local_fp.write(data)

    return local_path
