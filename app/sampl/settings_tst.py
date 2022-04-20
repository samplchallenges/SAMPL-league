# pylint: skip-file
from sampl.settings_dev import *  # lgtm [py/polluting-import]

CONTAINER_ENGINE = os.environ.get("CONTAINER_ENGINE", "docker")

DATABASES["default"]["NAME"] = "test_sampl"
