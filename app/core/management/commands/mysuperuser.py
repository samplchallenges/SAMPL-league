import os

from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand

USER_NAME = "admin"
EMAIL = "admin@example.com"
PASSWORD = os.environ["SUPER_USER_PASS"]


class Command(BaseCommand):
    def handle(self, *args, **options):  # pylint:disable=unused-argument
        User = get_user_model()
        # https://docs.djangoproject.com/en/3.2/topics/auth/customizing/#django.contrib.auth.get_user_model
        if not User.objects.filter(username=USER_NAME).exists():
            User.objects.create_superuser(USER_NAME, EMAIL, PASSWORD)
            self.stdout.write(self.style.SUCCESS("Superuser Created"))
        else:
            self.stdout.write(self.style.WARNING("Superuser exists"))
