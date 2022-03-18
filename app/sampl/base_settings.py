"""
Django settings for SAMPL project.
"""
import os
import subprocess

# Custom settings for SAMPL
VISUALIZE_DASK_GRAPH = False
DASK_SCHEDULER_URL = "localhost:8786"
ENABLE_REGISTRATION = False  # Disables self-registration
# End of custom settings for SAMPL

ALLOWED_HOSTS = []
# Application definition

DEFAULT_AUTO_FIELD = "django.db.models.AutoField"


INSTALLED_APPS = (
    "core.apps.CoreConfig",
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "crispy_forms",
)

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

ROOT_URLCONF = "sampl.urls"
LOGIN_REDIRECT_URL = "/"
LOGOUT_REDIRECT_URL = "/"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

WSGI_APPLICATION = "sampl.wsgi.application"


# Password validation
# https://docs.djangoproject.com/en/3.1/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]


# Internationalization
# https://docs.djangoproject.com/en/3.1/topics/i18n/

LANGUAGE_CODE = "en-us"

TIME_ZONE = "UTC"

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/3.1/howto/static-files/

STATIC_URL = "/static/"
STATIC_ROOT = "static"

CRISPY_TEMPLATE_PACK = "bootstrap4"

# container engine is either "docker" or "singularity"
CONTAINER_ENGINE = "docker"


ECR_BASE_URL = "103125031445.dkr.ecr.us-east-2.amazonaws.com"
ECR_SAMPLLEAGUE_URL = "103125031445.dkr.ecr.us-east-2.amazonaws.com/sampl-league"


def run_aws_login(container_engine):
    if container_engine == "docker":
        subprocess.run(
            [
                f"aws ecr get-login-password --region us-east-2 --profile sampl_pull | docker login --username AWS --password-stdin {ECR_BASE_URL}"
            ],
            shell=True,
            check=True,
        )
    elif container_engine == "singularity":
        os.environ["SINGULARITY_DOCKER_USERNAME"] = "AWS"
        aws_call = subprocess.run(
            ["aws ecr get-login-password --region us-east-2 --profile sampl_pull"],
            shell=True,
            capture_output=True,
            text=True,
            check=True,
        )
        os.environ["SINGULARITY_DOCKER_PASSWORD"] = aws_call.stdout
    else:
        raise Exception("Container Engine not yet implemented")


def run_aws_logout(container_engine):
    if container_engine == "singularity":
        os.environ["SINGULARITY_DOCKER_USERNAME"] = ""
        os.environ["SINGULARITY_DOCKER_PASSWORD"] = ""


LOGIN_TO_AWS = False
AWS_LOGIN_FUNCTION = run_aws_login
AWS_LOGOUT_FUNCTION = run_aws_logout
