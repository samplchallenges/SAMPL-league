from pathlib import Path
import typing

import docker  # type: ignore

from .utils import ContainerInstance, Engine, GUEST_OUTPUT_DIR

DOCKER_CONTAINER_TYPES = ["docker"]

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
        return self.container.status


class DockerEngine(Engine):
    _engine_name = "docker"

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
        if container_type not in DOCKER_CONTAINER_TYPES:
            raise ValueError(f"Container type: {container_type} not supported by Docker Engine")
        if aws_login_func:
            aws_login_func("docker")
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
