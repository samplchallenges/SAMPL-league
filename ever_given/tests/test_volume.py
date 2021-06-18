import os
import tempfile

import pytest

import ever_given.wrapper


def test_inputdir():
    test_mdlfile = "data/ChEBI_16716.mdl"
    test_data_dir = os.path.join(os.path.dirname(__file__), os.path.dirname(test_mdlfile))
    container_uri = "ghcr.io/robbason/calc-molwt:latest"
    command = f"--molfile /mnt/inputs/{os.path.basename(test_mdlfile)}"
    with tempfile.TemporaryDirectory() as tmpdir:
        inputdir = str(tmpdir)
        inputdir = test_data_dir
        result = ever_given.wrapper.run_container(container_uri, command, inputdir=inputdir)
    assert pytest.approx(float(result.strip()), 78.046950192)
