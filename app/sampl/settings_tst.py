# pylint: skip-file
import os

from sampl.settings_test import *  # lgtm [py/polluting-import]

CONTAINER_ENGINE = os.environ.get("CONTAINER_ENGINE", "docker")

DATABASES["default"]["NAME"] = "test_sampl"
