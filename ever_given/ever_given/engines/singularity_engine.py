from pathlib import Path

from spython.main import Client


from .utils import ContainerInstance, Engine, GUEST_OUTPUT_DIR


class SingularityContainerInstance(ContainerInstance):
    def __init__(self, container):
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


class SingularityEngine(Engine):
    _engine_name = "singularity"

    @classmethod
    def run_container(
        cls, container_uri, command_list, *, inputdir_map=None, output_dir=None
    ):
        bind_volumes = []
        for inputdir, guest_input_dir in inputdir_map.items():
            bind_volumes.append(f"{inputdir}:{guest_input_dir}:ro")
        if output_dir:
            output_dir = Path(output_dir).resolve()
            bind_volumes.append(f"{output_dir}:{GUEST_OUTPUT_DIR}:rw")
            command_list.extend(["--output-dir", GUEST_OUTPUT_DIR])

        sing_container = Client.run(
            image=container_uri,
            args=command_list,
            writable=True,
            contain=True,
            bind=bind_volumes,
            stream=True,
            # options, singularity_options,
            # nv=False, TODO: how to decide when to load Nvidia drivers?
            # network_disabled=True,
            # network_mode="none",
            remove=False,
            detach=True,
        )

        return SingularityContainerInstance(sing_container)
