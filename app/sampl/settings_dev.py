# pylint: skip-file
import os
from pathlib import Path

from .base_settings import *  # lgtm [py/polluting-import]

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent
MEDIA_ROOT = BASE_DIR / "media"

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True
# requires installing django-debug-toolbar
INSTALLED_APPS = (*INSTALLED_APPS, "debug_toolbar")
MIDDLEWARE = [
    "debug_toolbar.middleware.DebugToolbarMiddleware",
] + MIDDLEWARE
INTERNAL_IPS = ["127.0.0.1"]

# If in debug mode, have crispy-forms raise errors
# https://django-crispy-forms.readthedocs.io/en/latest/crispy_tag_forms.html#make-crispy-forms-fail-loud
CRISPY_FAIL_SILENTLY = not DEBUG

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = "7&gz49!+qmxwltqo*!_g$n+i)%qcn9c^%2kbwzlnwiuofj29!%"

if "RDS_DB_NAME" in os.environ:
    DATABASES = {
        "default": {
            "ENGINE": "django.db.backends.postgresql_psycopg2",
            "NAME": os.environ["RDS_DB_NAME"],
            "USER": os.environ["RDS_USERNAME"],
            "PASSWORD": os.environ["RDS_PASSWORD"],
            "HOST": os.environ["RDS_HOSTNAME"],
            "PORT": os.environ["RDS_PORT"],
        }
    }
else:
    DATABASES = {
        "default": {
            "ENGINE": "django.db.backends.sqlite3",
            "NAME": BASE_DIR / "db.sqlite3",
        }
    }

# LOGIN_TO_AWS = True
# AWS_LOGIN_FUNCTION = run_aws_login_singularity
