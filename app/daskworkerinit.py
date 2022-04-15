import django

django.setup()

from django.conf import settings

print(settings.CONTAINER_ENGINE)
