import logging
from pathlib import Path
import typing
import subprocess

import docker  # type: ignore

from .utils import ContainerInstance, Engine, GUEST_OUTPUT_DIR

logger = logging.getLogger(__name__)


class DockerContainerInstance(ContainerInstance):
    def __init__(self, client, container):
        self.client = client
        self.container = container

    def logs(self, *, want_stdout, want_stderr):
        # Docker gives us bytes; convert to text with UTF-8 encoding
        for row in self.container.logs(
            stdout=want_stdout, stderr=want_stderr, stream=True
        ):
            yield row.decode()

    def reload(self):
        self.container.reload()

    def kill(self):
        self.container.kill()

    def remove(self):
        self.container.remove()

    def status(self):
        if self.container.status == "exited":
            if self.container.attrs["State"]["ExitCode"] == 0:
                return ContainerInstance.SUCCESS
            
            return ContainerInstance.FAILURE
        if self.container.status == "running":
            return ContainerInstance.RUNNING
        raise ValueError(f"Unknown docker container status: {self.container.status}")


class DockerEngine(Engine):
    _engine_name = "docker"
    _valid_container_types = ["docker"]

    @classmethod
    def run_container(
        cls,
        container_type: str,
        container_uri: str,
        command_list: typing.List[str],
        *,
        inputdir_map: typing.Dict[str, str] = None,
        output_dir: str = None,
        aws_login_func=None,
    ) -> DockerContainerInstance:
        cls.validate_common_arguments(container_type, aws_login_func)
        command = " ".join(command_list)
        client = docker.from_env()
        volumes = {}
        if inputdir_map is not None:
            for inputdir, guest_input_dir in inputdir_map.items():
                volumes[str(inputdir)] = {"bind": str(guest_input_dir), "mode": "ro"}
        if output_dir:
            output_dir = str(Path(output_dir).resolve())
            volumes[output_dir] = {"bind": str(GUEST_OUTPUT_DIR), "mode": "rw"}
            command = f" --output-dir {GUEST_OUTPUT_DIR} {command}"
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Docker command: %s", command)
            logger.debug("Volumes: %s", volumes)
        docker_container = client.containers.run(
            container_uri,
            command,
            volumes=volumes,
            network_disabled=True,
            network_mode="none",
            remove=False,
            detach=True,
        )
        return DockerContainerInstance(client, docker_container)

    @classmethod
    def pull_container(cls, container_uri, container_type, save_path=None, aws_login_func=None):
        # pylint: disable=unused-argument
        cls.validate_common_arguments(container_type, aws_login_func)

        pull_cmd = ['docker', 'pull', container_uri]

        ended_proc = subprocess.run(pull_cmd, capture_output=True)
        code = ended_proc.returncode
        stdout = ended_proc.stdout.decode("utf-8")
        stderr = ended_proc.stderr.decode("utf-8")
        return code == 0, stdout, stderr
        
