"""
Shared classes and constants used by docker and singularity
"""
import abc
from pathlib import Path
import typing

GUEST_OUTPUT_DIR = Path("/mnt") / "outputs"


class ContainerInstance(abc.ABC):
    SUCCESS = "success"
    FAILURE = "failure"
    RUNNING = "running"

    @abc.abstractmethod
    def logs(
        self, *, want_stdout: bool, want_stderr: bool
    ) -> typing.Generator[str, None, None]:
        pass

    @abc.abstractmethod
    def reload(self) -> None:
        pass

    @abc.abstractmethod
    def kill(self) -> None:
        pass

    @abc.abstractmethod
    def remove(self) -> None:
        pass

    @abc.abstractmethod
    def status(self) -> str:
        pass


class Engine(abc.ABC):
    """
    Base class used by docker and singularity enngines
    """

    _engine_name: str
    _valid_container_types: typing.List[str]

    @classmethod
    def name(cls) -> str:
        if cls._engine_name is None:
            raise ValueError("Engine name is not set")
        return cls._engine_name

    @classmethod
    @abc.abstractmethod
    def run_container(
        cls,
        container_type: str,
        container_uri: str,
        command_list: typing.List[str],
        *,
        inputdir_map: typing.Dict[str, str] = None,
        output_dir: str = None,
        aws_login_func: typing.Callable[[str], None]=None
    ) -> ContainerInstance:
        pass

    @classmethod
    @abc.abstractmethod
    def pull_container(cls, container_uri: str, container_type: str, save_path: str=None, aws_login_func: typing.Callable[[str], None]=None) -> typing.Tuple[bool, str, str]:
        """
        Returns False, stdout, stderr on failure, True, stdout, stderr on success.
        Note this is an inversion of UNIX, where 0 exit code is success (which is falsy in Python)"""

    @classmethod
    def validate_common_arguments(cls, container_type: str, aws_login_func: typing.Callable[[str], None]) -> None:
        if container_type not in cls._valid_container_types:
            raise ValueError(f"Container type: {container_type} not supported by {cls._engine_name} Engine")
        if aws_login_func:
            aws_login_func(cls._engine_name)
