import os

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "console2": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
        },
        "console": {
            "level": "DEBUG",
            "class": "logging.FileHandler",
            "filename": "test_logging.log",
        },
    },
    "loggers": {
        "django": {
            "handlers": ["console"],
            "level": os.getenv("DJANGO_LOG_LEVEL", "DEBUG"),
            "propagate": False,
        },
        "core": {
            "handlers": ["console"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "DEBUG"),
            "propagate": False,
        },
        "referee": {
            "handlers": ["console"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "DEBUG"),
            "propagate": False,
        },
        "ever_given": {
            "handlers": ["console"],
            "level": os.getenv("SAMPL_LOG_LEVEL", "DEBUG"),
            "propagate": False,
        },
    },
}
