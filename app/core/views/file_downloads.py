from django.contrib.auth.decorators import login_required
from django import http
from django.conf import settings
from django.shortcuts import get_object_or_404

from .. import filecache
from .. import models
#import ContainerArg, InputValue


def _respond_file(value):
    if settings.DEFAULT_FILE_STORAGE == "django.core.files.storage.FileSystemStorage":
        local_path = filecache.ensure_local_copy(value)
        return http.FileResponse(open(local_path, "rb"), as_attachment=True)

    return http.HttpResponseRedirect(value.url)


@login_required
def download_output_file(request, pk):  # pylint: disable=unused-argument
    prediction = get_object_or_404(
        models.Prediction,
        evaluation__submission_run__submission__user=request.user,
        pk=pk)

    return _respond_file(prediction.value)


@login_required
def download_input_file(request, pk):  #  pylint: disable=unused-argument
    input_value = get_object_or_404(models.InputValue,
        input_element__is_public=True,
        value_type__content_type__model="filevalue",
        pk=pk,
    )

    return _respond_file(input_value.value)


@login_required
def download_container_arg_file(request, pk):
    container_arg = get_object_or_404(models.ContainerArg,
        container__user=request.user, pk=pk)

    return _respond_file(container_arg.file_value)
