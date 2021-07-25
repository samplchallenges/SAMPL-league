"""
When we have files stored on S3, we want to seamlessly manage local copies for submission to containers.
"""
import logging
import os.path

import boto3

from django.conf import settings


BUFSIZE = 10*1024*1024  # Break big files into 10MB chunks
# BUFSIZE = 100

logger = logging.getLogger(__name__)


def ensure_local_copy(field_file):
    """
    Should be a no-op if using default Django storage
    """
    local_path = os.path.join(settings.MEDIA_ROOT, field_file.name)
    if not os.path.exists(local_path):
        logger.debug("Creating local cache %s", local_path)
        os.makedirs(os.path.dirname(local_path))
        with open(local_path, 'wb') as local_fp:
            while True:
                data = field_file.read(BUFSIZE)
                if not data:
                    break
                local_fp.write(data)
                logger.debug("Writing chunk")

    return local_path
