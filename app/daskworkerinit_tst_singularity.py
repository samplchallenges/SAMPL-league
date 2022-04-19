import os
import django

os.environ["DJANGO_SETTINGS_MODULE"] = "sampl.settings_tst_singularity"

django.setup()
