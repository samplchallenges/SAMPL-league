from asyncore import poll
import io
from pathlib import Path
import subprocess
import typing
import shlex
import os

from .utils import ContainerInstance, Engine, GUEST_OUTPUT_DIR

AFTER_TERMINATE_WAIT = 5  # Seconds to wait after 'terminate' signal before 'kill'


class SingularityContainerInstance(ContainerInstance):
    def __init__(self, process: subprocess.Popen, container_uri):
        self.process = process
        self.container_uri = container_uri

    def logs(self, *, want_stdout: bool, want_stderr: bool):
        if want_stdout and want_stderr:
            raise ValueError("Can't have want_stdout and want_stderr both true")
        pipe: io.BufferedReader
        if want_stdout:
            pipe = typing.cast(io.BufferedReader, self.process.stdout)
        elif want_stderr:
            pipe = typing.cast(io.BufferedReader, self.process.stderr)
        else:
            raise ValueError("Can't have want_stdout and want_stderr both false")
        
        for line in iter(pipe.readline, b''):
            yield line.decode('utf-8')

    def reload(self):
        pass  # don't need this for status, right?

    def kill(self):
        self.process.terminate()
        try:
            self.process.wait(AFTER_TERMINATE_WAIT)
        except subprocess.TimeoutExpired:
            self.process.kill()

    def remove(self):
        pass
        # raise Exception("Not yet implemented")
        # self.container.remove()

        # container_sif = os.path.basename(self.container_uri)
        # command = f"singularity cache clean "
        # this removes the entire cache, singularity 3.7.2 does not have a --name option
        # to remove specific containers

    def status(self):
        poll_status = self.process.poll()
        if poll_status is None:
            return ContainerInstance.RUNNING
        
        if poll_status == 0:
            return ContainerInstance.SUCCESS
        
        return ContainerInstance.FAILURE


class SingularityEngine(Engine):
    _engine_name = "singularity"
    uri_prefixes: typing.Dict[str, str] = {
        "docker": "docker://",
        "singularity_remote": "shub://",
        "singularity_local": ""}
    _valid_container_types = list(uri_prefixes.keys())

    @classmethod
    def run_container(
        cls,
        container_type,
        container_uri,
        command_list,
        *,
        inputdir_map=None,
        output_dir=None,
        aws_login_func=None,
    ):
        cls.validate_common_arguments(container_type, aws_login_func)

        bind_volumes = []
        for inputdir, guest_input_dir in inputdir_map.items():
            bind_volumes.append(f"{inputdir}:{guest_input_dir}")  #:ro")
        if output_dir:
            output_dir = Path(output_dir).resolve()
            bind_volumes.append(f"{output_dir}:{GUEST_OUTPUT_DIR}")  #:rw")
            command_list.extend(["--output-dir", GUEST_OUTPUT_DIR])

        command = _build_singularity_command(bind_volumes, container_uri, command_list)
        process = subprocess.Popen(
            command, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return SingularityContainerInstance(process, container_uri)

    @classmethod
    def make_uri(cls, container_uri, container_type):

        try:
            return f"{cls.uri_prefixes[container_type]}{container_uri}"
        except KeyError:
            raise Exception(
                f"Container Type {container_type} not implemented,"
                f" valid options are {cls.uri_prefixes.keys()}")


    @classmethod
    def pull_container(cls, container_uri, container_type, save_path=None, aws_login_func=None):
        cls.validate_common_arguments(container_type, aws_login_func)

        uri = cls.make_uri(container_uri, container_type)
        if save_path:   
            pull_cmd = ["singularity", "pull", "-F", save_path, uri]
        else:
            pull_cmd = ["singularity", "pull", "-F", uri]
        ended_proc = subprocess.run(pull_cmd, capture_output=True)
        code = ended_proc.returncode
        stdout = ended_proc.stdout.decode("utf-8")
        stderr = ended_proc.stderr.decode("utf-8")
        return code == 0, stdout, stderr
            

def _build_singularity_command(bind_volumes, container_uri, command_list):
    #  TODO: how to decide when to load Nvidia drivers?
    singularity_cmd = [
        "singularity",
        "run",
    ]
    if bind_volumes and len(bind_volumes) >= 1:
        singularity_cmd.append("--bind")
        bind_str = ""
        for bind in bind_volumes:
            bind_str += bind + ","
        singularity_cmd.append(bind_str[:-1])

    singularity_cmd.extend([container_uri] + command_list)

    return singularity_cmd
