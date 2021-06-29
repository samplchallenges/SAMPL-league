import copy
from pathlib import Path
import os
import os.path
import shlex
import tempfile

import docker

GUEST_INPUT_DIR = Path("/mnt") / "inputs"
GUEST_OUTPUT_DIR = Path("/mnt") / "outputs"

def _guest_input_path(filename):
    return GUEST_INPUT_DIR / Path(filename).name


def _prepare_commandline(command, args_dict):
    return command + " ".join(
        [f"--{key} {shlex.quote(value)}" for key, value in args_dict.items()]
    )


def _parse_output(raw_text):
    result = {}
    for line in raw_text.decode("utf-8").splitlines():
        lineparts = line.split(maxsplit=1)
        if len(lineparts) == 2:
            key, value = lineparts
            result[key] = value
        else:
            raise ValueError(f"Cannot parse output {line}, needs KEY VALUE format")
    return result


def _copyfile(src, dest):
    # TODO: better way to handle big files
    # can we avoid copying entirely? Mount multiple dirs read-only
    BUFSIZE = 10*1024*2
    while True:
        buf = src.read(BUFSIZE)
        if buf:
            dest.write(buf)
        else:
            break


def run(container_uri, command="", *, has_output_files,
        file_kwargs, kwargs):
    """
    kwargs will be passed to container as --key=value
    where value will be shell escaped
    file_kwargs must be a dict of key to Python file-like object
    the underlying files will be copied to a temporary directory,
    which will be made available to the container as /mnt/inputs (readonly)
    if has_output_files is True, a read-write directory will be passed
    to the container
    command is optional
    """
    final_kwargs = copy.deepcopy(kwargs)
    with tempfile.TemporaryDirectory() as tmpdir:
        dirpath = Path(str(tmpdir))
        inputdir = dirpath / "inputs"
        inputdir.mkdir()
        if has_output_files:
            outputdir = dirpath / "output"
            outputdir.mkdir()
            scoringdir = dirpath / "scoring"
            scoringdir.mkdir()
            (scoringdir / "predictions").mkdir()
            (scoringdir / "answers").mkdir()
        else:
            outputdir = None
            scoringdir = None


        for key, fp in file_kwargs.items():
            final_kwargs[key] = _guest_input_path(fp.name)
            #with open(inputdir / fp.name, "wb") as local_fp:
            #    _copyfile(fp, local_fp)
            #shutil.copy(fp.path, inputdir)

        final_command = _prepare_commandline(command, final_kwargs)

        result = run_container(
            container_uri, final_command, inputdir=inputdir, outputdir=outputdir
        )
        for item in _parse_output(result).items():
            yield item



def run_container(container_uri, command, inputdir=None, outputdir=None):
    client = docker.from_env()
    volumes = {}
    if inputdir:
        volumes[inputdir] = {
            'bind': str(GUEST_INPUT_DIR),
            'mode': 'ro'}
    if outputdir:
        volumes[outputdir] = {
            'bind': str(GUEST_OUTPUT_DIR),
            'mode': 'rw'}
        command = f" --output-dir {GUEST_OUTPUT_DIR} {command}"

    result = client.containers.run(
        container_uri,
        command,
        volumes=volumes,
        network_disabled=True,
        network_mode="none",
        remove=True,
    )
    return result
