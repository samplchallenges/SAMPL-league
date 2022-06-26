from django.apps import AppConfig


class CoreConfig(AppConfig):
    name = "core"

    def ready(self):
        from . import signals  # pylint: disable=unused-import

        # NB: this may cause signals to run multiple times in testing
        # https://docs.djangoproject.com/en/4.0/topics/signals/#preventing-duplicate-signals
