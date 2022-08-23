LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "console": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
        },
    },
    "loggers": {
        "django": {
            "handlers": ["console"],
            "level": "DEBUG",
            "propagate": False,
        },
        "core": {
            "handlers": ["console"],
            "level": "DEBUG",
            "propagate": False,
        },
        "referee": {
            "handlers": ["console"],
            "level": "DEBUG",
            "propagate": False,
        },
        "ever_given": {
            "handlers": ["console"],
            "level": "DEBUG",
            "propagate": False,
        },
    },
}
