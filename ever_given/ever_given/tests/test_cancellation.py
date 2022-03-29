import time
import os
import pytest
import subprocess
import shlex

import ever_given.wrapper
from ever_given.log_processing import QUEUE_WAIT_SECONDS, CancelledException

@pytest.mark.parametrize(["container_engine"], [["docker"], ["singularity"]])
def test_cancellation(container_engine):
    # Use a slow container so we have time to cancel it
    container_uri = "ghcr.io/megosato/logging-example:latest"
    container_type = "docker"
    if container_engine == "singularity":
        command = ["singularity", "pull", f"docker://{container_uri}"]
        subprocess.run(command, check=True)
        container_uri = "logging-example_latest.sif"
        container_type = "singularity_local"
        assert os.path.exists(container_uri)
    if container_engine == "docker":
        command = ["docker", "pull", container_uri]
        subprocess.run(command, shell=True, check=True)
    kwargs = {"smiles": "c1cccnc1"}

    start_at = time.time()

    with pytest.raises(CancelledException):
        results = {
            key: value
            for key, value in ever_given.wrapper.run(
                container_uri,
                kwargs=kwargs,
                file_kwargs={},
                container_type=container_type,
                engine_name=container_engine,
                cancel_requested_func=lambda: True,
            )
        }
    end_at = time.time()
    # Should cancel in not much more time than QUEUE_WAIT_SECONDS * 2
    # in log_processing
    assert end_at - start_at < QUEUE_WAIT_SECONDS * 2 + 2



@pytest.mark.parametrize(["container_engine"], [["docker"], ["singularity"]])
def test_no_cancellation(container_engine):
    container_uri = "ghcr.io/megosato/calc-molwt:latest"
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri,
            kwargs=kwargs,
            file_kwargs={},
            container_type="docker",
            engine_name=container_engine,
            cancel_requested_func=lambda: False,
        )
    }
