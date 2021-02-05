import os
from django.core.management.base import BaseCommand
from django.contrib.auth.models import User

USER_NAME = "admin"
EMAIL = "admin@example.com"
PASSWORD = os.environ["SUPER_USER_PASS"]

class Command(BaseCommand):
    def handle(self, *args, **options):
        if not User.objects.filter(username=USER_NAME).exists():
            User.objects.create_superuser(USER_NAME,
                                          EMAIL,
                                          PASSWORD)
