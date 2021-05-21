import sys
import os

import django

from django.conf import settings

# from sampl import settings_tst

os.environ["DJANGO_SETTINGS_MODULE"] = "sampl.settings_tst"

# dir_path = os.path.dirname(os.path.realpath(__file__))
# sys.path.insert(0, dir_path)
if not settings.configured:
    #    settings.configure(default_settings=settings_tst)
    django.setup()

print(settings.DATABASES)
