from pathlib import Path

import docker


from .utils import ContainerInstance, Engine, GUEST_OUTPUT_DIR


class DockerContainerInstance(ContainerInstance):
    def __init__(self, client, container):
        self.client = client
        self.container = container

    def logs(self, *, stdout, stderr):
        return self.container.logs(stdout=stdout, stderr=stderr, stream=True)

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
        cls, container_uri, command_list, *, inputdir_map=None, output_dir=None
    ):
        command = " ".join(command_list)
        client = docker.from_env()
        volumes = {}
        for inputdir, guest_input_dir in inputdir_map.items():
            volumes[str(inputdir)] = {"bind": str(guest_input_dir), "mode": "ro"}
        if output_dir:
            output_dir = Path(output_dir).resolve()
            volumes[str(output_dir)] = {"bind": str(GUEST_OUTPUT_DIR), "mode": "rw"}
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
