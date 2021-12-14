"""
Shared classes and constants used by docker and singularity
"""
import abc
from pathlib import Path


GUEST_OUTPUT_DIR = Path("/mnt") / "outputs"


class Engine(abc.ABC):
    """
    Base class used by docker and singularity enngines
    """

    _engine_name = None

    @classmethod
    @property
    def name(cls):
        return cls._engine_name

    @classmethod
    @abc.abstractmethod
    def run_container(
        cls, container_uri, command, *, inputdir_map=None, output_dir=None
    ):
        pass


class ContainerInstance(abc.ABC):
    @abc.abstractmethod
    def logs(self, container):
        pass

    @abc.abstractmethod
    def reload(self):
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
