from sampl.settings_dev import *  # lgtm [py/polluting-import]

DATABASES["default"] = {
    "ENGINE": "django.db.backends.postgresql_psycopg2",
    "NAME": "test_sampl",
    "USER": "sampl",
    "HOST": "localhost"
}
