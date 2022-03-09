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
        string = ""
        for ln in pipe.readlines():
            string += ln.decode('utf-8')
        return string

    def reload(self):
        pass  # don't need this for status, right?

    def kill(self):
        self.process.terminate()
        try:
            self.process.wait(AFTER_TERMINATE_WAIT)
        except subprocess.TimeoutExpired:
            self.process.kill()

    def remove(self):
        print(self.container_uri)
        
        #raise Exception("Not yet implemented")
        # self.container.remove()

        #container_sif = os.path.basename(self.container_uri)
        #command = f"singularity cache clean -N {container_sif}"

        #output = subprocess.run(command, shell=True, capture_output=True)
        #print(output)
         
    def status(self):
        #print("status:", self.process.poll())
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
            bind_volumes.append(f"{inputdir}:{guest_input_dir}") #:ro")
        if output_dir:
            output_dir = Path(output_dir).resolve()
            bind_volumes.append(f"{output_dir}:{GUEST_OUTPUT_DIR}") #:rw")
            command_list.extend(["--output-dir", GUEST_OUTPUT_DIR])


        command = _build_singularity_command(bind_volumes, container_uri, command_list)
        print(command)
        process = subprocess.Popen(
            shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return SingularityContainerInstance(process, container_uri)



def _build_singularity_command(bind_volumes, container_uri, command_list):
    #  TODO: how to decide when to load Nvidia drivers?
    singularity_cmd = "singularity run --bind "
    for bind in bind_volumes:
        singularity_cmd += bind + ","
    singularity_cmd = singularity_cmd[:-1] + " "
    singularity_cmd += container_uri + " "
    for command in command_list:
        singularity_cmd += str(command) + " "
    return singularity_cmd
    #raise Exception("not implemented yet")

