import copy
import logging
import sys
from pathlib import Path
import typing

from .engines import GUEST_OUTPUT_DIR, REGISTERED_ENGINES

from . import log_processing
from .utils import LogHandlerBase


GUEST_INPUT_DIR = Path("/mnt") / "inputs"
logger = logging.getLogger(__name__)


class ContainerFailedException(Exception):
    "If the container process exits with a non-zero return value"


def _guest_input_path(filename: str) -> Path:
    return GUEST_INPUT_DIR / Path(filename).name


def prepare_command_list(command: str, args_dict: typing.Dict[str, str]):
    command_list = [command] if command else []
    for key, value in args_dict.items():
        command_list.extend([f"--{key}", str(value)]) 
    return command_list


def _parse_output(
    host_output_path: typing.Optional[str],
    raw_text: str,
    output_file_keys: typing.Optional[typing.List[str]],
) -> typing.Dict[str, str]:
    result = {}
    for line in raw_text.splitlines():
        lineparts = line.split(maxsplit=1)
        if len(lineparts) == 2:
            key, value = lineparts
            value = value.strip()
            if output_file_keys and key in output_file_keys:
                if host_output_path is None:
                    raise ValueError("output path not set but output files expected")
                print(host_output_path, key, value)
                result[key] = _convert_guest_to_host_path(host_output_path, value)
            else:
                result[key] = value

    return result


def _convert_guest_to_host_path(host_path: str, filepath: str) -> str:
    abs_filepath = GUEST_OUTPUT_DIR / filepath
    relative_path = abs_filepath.relative_to(GUEST_OUTPUT_DIR)
    # so we don't care whether container outputs relative or absolute paths
    return str(Path(host_path) / relative_path)


def _convert_file_kwargs(
    file_kwargs: typing.Dict[str, str]
) -> typing.Tuple[typing.Dict[Path, Path], typing.Dict[str, str]]:
    """
    file_kwargs has key: path, where path is on the host
    generate a dict of host: guest path mappings that can be used to generate bind commands and return a file_kwargs pointing to the files with corresponding paths that are accessible to the container
    """
    dirpaths = {}
    final_file_kwargs = {}
    for idx, (key, pathname) in enumerate(file_kwargs.items()):
        path = Path(pathname).resolve()
        dirpath = path.parent
        basename = path.name
        if dirpath not in dirpaths:
            dirpaths[dirpath] = Path("/mnt") / f"inputs{idx}"

        key = str(key)
        final_file_kwargs[key] = str(dirpaths[dirpath] / basename)
    return dirpaths, final_file_kwargs


def _get_container_uri(container_uri, container_type, engine_name):
    if engine_name == "singularity":
        return REGISTERED_ENGINES["singularity"].make_uri(container_uri, container_type)
    else:
        return container_uri


def pull_container(container_uri, container_type, engine_name, save_path, aws_login_func=None):
    return REGISTERED_ENGINES[engine_name].pull_container(container_uri, container_type, save_path, aws_login_func)


def run(
    container_uri: str,
    command: str = "",
    *,
    file_kwargs: typing.Dict[str, str],
    kwargs: typing.Dict[str, str],
    container_type: str = "docker",
    engine_name: str = "docker",
    output_dir: str = None,
    output_file_keys: typing.List[str] = None,
    log_handler: LogHandlerBase = None,
    cancel_requested_func: typing.Callable[[], bool] = None,
    aws_login_func=None,
):
    """
    kwargs will be passed to container as --key=value
    where value will be shell escaped
    file_kwargs must be a dict of key to Python file-like object
    the underlying files will be made available to the container as bound mount points (readonly).
    Note that we must ensure the directories thus made available to the container only have challenge-wide data, no outputs from other submissions!
    output_dir, if not None, will be mounted r/w for container to write outputs
    output_file_keys is a set of keys for file output types. Their values will be mapped from paths on the guest to paths on the host.
    command is optional
    Iterate over key/values of results.
    """
    if log_handler is None:
        log_handler = LogHandlerBase()

    inputdir_map, final_file_kwargs = _convert_file_kwargs(file_kwargs)
    final_kwargs = copy.deepcopy(kwargs)
    final_kwargs.update(final_file_kwargs)
    command_list = prepare_command_list(command, final_kwargs)
    logger.debug("Command: %s", command_list)
    container_uri = _get_container_uri(container_uri, container_type, engine_name)
    logger.debug("Container URI: %s", container_uri)

    running_container = REGISTERED_ENGINES[engine_name].run_container(
        container_type,
        container_uri,
        command_list,
        inputdir_map=inputdir_map,
        output_dir=output_dir,
        aws_login_func=aws_login_func,
    )

    try:
        result = log_processing.process_messages(
            running_container, log_handler, cancel_requested_func
        )

        yield from _parse_output(output_dir, result, output_file_keys).items()

        running_container.reload()
        status = running_container.status()
        print("container status is ", status)
        if status == running_container.RUNNING:
            log_handler.handle_stderr("Killing running container\n")
            running_container.kill()
            running_container.reload()
            status = running_container.status()
            log_handler.handle_stderr(
                f"After killing running container, status is {status}\n"
            )
            print(f"Container status is {status}", file=sys.stderr)
        if status == running_container.FAILURE:
            # if log_handler.cls:
            #     bev = log_handler.cls.objects.get(id=log_handler.instance_id)
            #     print(container_uri)
            #     print(bev.log_stderr)
            # for debugging, above four lines can be uncommented
            raise ContainerFailedException()
    finally:
        running_container.remove()
