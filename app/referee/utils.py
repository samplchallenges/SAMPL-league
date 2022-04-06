from django.conf import settings


def get_aws_credential_function(container_uri):
    if container_uri.startswith(settings.ECR_BASE_URL):
        return settings.AWS_LOGIN_FUNCTION
    else:
        return settings.AWS_LOGOUT_FUNCTION
