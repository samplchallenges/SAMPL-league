import os

from .base_settings import *  # lgtm [py/polluting-import]

DEBUG = False
SECRET_KEY = os.environ["SAMPL_SECRET_KEY"]
ALLOWED_HOSTS = [
    "app.samplchallenges.org",
    "samplmvp-env.eba-rhcwa63p.us-east-2.elasticbeanstalk.com",
]

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
