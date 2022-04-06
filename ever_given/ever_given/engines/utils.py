"""
Shared classes and constants used by docker and singularity
"""
import abc
from pathlib import Path
import typing

GUEST_OUTPUT_DIR = Path("/mnt") / "outputs"


class ContainerInstance(abc.ABC):
    @abc.abstractmethod
    def logs(
        self, *, want_stdout: bool, want_stderr: bool
    ) -> typing.Generator[str, None, None]:
        pass

    @abc.abstractmethod
    def reload(self):
        pass

    @abc.abstractmethod
    def kill(self):
        pass

    @abc.abstractmethod
    def remove(self):
        pass

    @abc.abstractmethod
    def status(self):
        pass


class Engine(abc.ABC):
    """
    Base class used by docker and singularity enngines
    """

    _engine_name: typing.Optional[str] = None

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
        aws_login_func=None
    ) -> ContainerInstance:
        pass
