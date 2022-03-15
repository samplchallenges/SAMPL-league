import io
from pathlib import Path
import subprocess
import typing

from .utils import ContainerInstance, Engine, GUEST_OUTPUT_DIR

AFTER_TERMINATE_WAIT = 5  # Seconds to wait after 'terminate' signal before 'kill'


class SingularityContainerInstance(ContainerInstance):
    def __init__(self, process: subprocess.Popen):
        self.process = process

    def logs(self, *, want_stdout: bool, want_stderr: bool):
        if want_stdout and want_stderr:
            raise ValueError("Can't have want_stdout and want_stderr both true")
        pipe: io.TextIOBase
        if want_stdout:
            pipe = typing.cast(io.TextIOBase, self.process.stdout)
        elif want_stderr:
            pipe = typing.cast(io.TextIOBase, self.process.stderr)
        else:
            raise ValueError("Can't have want_stdout and want_stderr both false")
        return pipe.readline()

    def reload(self):
        pass  # don't need this for status, right?

    def kill(self):
        self.process.terminate()
        try:
            self.process.wait(AFTER_TERMINATE_WAIT)
        except subprocess.TimeoutExpired:
            self.process.kill()

    def remove(self):
        raise Exception("Not yet implemented")
        # self.container.remove()

    def status(self):
        return self.process.poll()


class SingularityEngine(Engine):
    _engine_name = "singularity"

    @classmethod
    def run_container(
        cls,
        container_uri,
        command_list,
        *,
        inputdir_map=None,
        output_dir=None,
        aws_login_func=None,
    ):
        if aws_login_func:
            aws_login_func()
        bind_volumes = []
        for inputdir, guest_input_dir in inputdir_map.items():
            bind_volumes.append(f"{inputdir}:{guest_input_dir}:ro")
        if output_dir:
            output_dir = Path(output_dir).resolve()
            bind_volumes.append(f"{output_dir}:{GUEST_OUTPUT_DIR}:rw")
            command_list.extend(["--output-dir", GUEST_OUTPUT_DIR])

        command = _build_singularity_command()
        process = subprocess.Popen(
            command, stdin=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return SingularityContainerInstance(process)


def _build_singularity_command():
    #  TODO: how to decide when to load Nvidia drivers?
    raise Exception("not implemented yet")
