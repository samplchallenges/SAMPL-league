from django.contrib.auth.decorators import login_required
from django.http import FileResponse

from .. import filecache
from ..models import InputValue, SubmissionArg


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
def download_submission_arg_file(request, pk):
    submission_arg = SubmissionArg.objects.get(submission__user=request.user, pk=pk)

    return _respond_local_file(submission_arg.file_value)
