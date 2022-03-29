import io
from pathlib import Path
import subprocess
import typing
import shlex
import os

from .utils import ContainerInstance, Engine, GUEST_OUTPUT_DIR

AFTER_TERMINATE_WAIT = 5  # Seconds to wait after 'terminate' signal before 'kill'
SINGULARITY_CONTAINER_TYPES = ["docker", "singularity_remote", "singularity_local"]

class SingularityContainerInstance(ContainerInstance):
    def __init__(self, process: subprocess.Popen, container_uri):
        self.process = process
        self.container_uri = container_uri

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

        # FIXED: return pipe.readline() 
        #   * returns an integer that needed to be decoded .decode('utf-8')
        #   * only returned the first line
        log_string_list = typing.cast(typing.List[str], [])

        while True:
            line = typing.cast(bytes, pipe.readline())
            if not line:
                break
            log_string_list.append(line.decode('utf-8'))
        return log_string_list

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
        #raise Exception("Not yet implemented")
        # self.container.remove()

        #container_sif = os.path.basename(self.container_uri)
        #command = f"singularity cache clean "
        #this removes the entire cache, singularity 3.7.2 does not have a --name option
        # to remove specific containers
         
    def status(self):
        return self.process.poll()


class SingularityEngine(Engine):
    _engine_name = "singularity"

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
        # Check the container type
        if container_type not in SINGULARITY_CONTAINER_TYPES:
            raise ValueError(f"Container type {container_type} not supported by singularity engine")
        if aws_login_func:
            aws_login_func("singularity")
        bind_volumes = []
        for inputdir, guest_input_dir in inputdir_map.items():
            bind_volumes.append(f"{inputdir}:{guest_input_dir}") #:ro")
        if output_dir:
            output_dir = Path(output_dir).resolve()
            bind_volumes.append(f"{output_dir}:{GUEST_OUTPUT_DIR}") #:rw")
            command_list.extend(["--output-dir", GUEST_OUTPUT_DIR])


        command = _build_singularity_command(bind_volumes, container_uri, command_list)
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return SingularityContainerInstance(process, container_uri)

    @classmethod
    def make_uri(cls, container_uri, container_type):
        if container_type == "docker":
            return f"docker://{container_uri}"
        elif container_type == "singularity_remote":
            return f"shub://{container_uri}"
        elif container_type == "singularity_local":
            return container_uri
        else:
            raise Exception("Container Type Not Implemented")

def _build_singularity_command(bind_volumes, container_uri, command_list):
    #  TODO: how to decide when to load Nvidia drivers?
    singularity_cmd = ["singularity", "run",]
    if bind_volumes and len(bind_volumes) >= 1:
        singularity_cmd.append("--bind")
        bind_str = ""
        for bind in bind_volumes:
            bind_str += bind + ","
        singularity_cmd.append(bind_str[:-1])
    
    singularity_cmd.extend([container_uri] + command_list)
    
    return singularity_cmd

