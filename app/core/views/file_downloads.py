from django.contrib.auth.decorators import login_required
from django.http import FileResponse

from .. import filecache
from ..models import ContainerArg, InputValue


def _respond_local_file(value):
    local_path = filecache.ensure_local_copy(value)
    response = FileResponse(open(local_path, "rb"))
    return response


@login_required
def download_input_file(request, pk):  #  pylint: disable=unused-argument
    input_value = InputValue.objects.get(
        input_element__is_public=True,
        value_type__content_type__model="filevalue",
        pk=pk,
    )

    return _respond_local_file(input_value.value)


@login_required
def download_container_arg_file(request, pk):
    container_arg = ContainerArg.objects.get(container__user=request.user, pk=pk)

    return _respond_local_file(container_arg.file_value)
