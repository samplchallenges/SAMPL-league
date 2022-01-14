"""
WSGI config for SAMPL project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/3.1/howto/deployment/wsgi/
"""

import json
import os

from django.core.wsgi import get_wsgi_application

with open("/etc/sampl_rds_config.json") as fp:
    rds_config = json.load(fp)

for key, value in rds_config.items():
    os.environ.setdefault(key, value)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "sampl.settings_prod")

application = get_wsgi_application()
