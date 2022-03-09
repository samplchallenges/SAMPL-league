# singularity tests
import os.path
import tempfile
import time

import pytest

import ever_given.wrapper
from ever_given.log_processing import QUEUE_WAIT_SECONDS, CancelledException

def test_run_inputfile_only():
    test_mdlfile_rel = "data/ChEBI_16716.mdl"
    test_mdlfile_abs = os.path.join(os.path.dirname(__file__), test_mdlfile_rel)
    cached_container = os.path.join(os.path.dirname(__file__), "data/calc-molwt_latest.sif")
    repo_container = "docker://ghcr.io/megosato/calc-molwt"
    container_uri = cached_container if os.path.exists(cached_container) else repo_container
    file_kwargs = {"molfile": test_mdlfile_abs}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri, 
            engine_name="singularity", 
            kwargs={}, 
            file_kwargs=file_kwargs
        )
    }
    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(78.046950192)

def test_cancellation():
    # Use a slow container so we have time to cancel it
    cached_container = os.path.join(os.path.dirname(__file__), "data/logging-example_latest.sif")
    repo_container = "docker://ghcr.io/megosato/logging-example"
    container_uri = cached_container if os.path.exists(cached_container) else repo_container

    kwargs = {"smiles": "c1cccnc1"}
    start_at = time.time()
    with pytest.raises(CancelledException):
        results = {
            key: value
            for key, value in ever_given.wrapper.run(
                container_uri,
                engine_name="singularity",
                kwargs=kwargs,
                file_kwargs={},
                cancel_requested_func=lambda: True,
            )
        }
    end_at = time.time()
    # Should cancel in not much more time than QUEUE_WAIT_SECONDS * 2
    # in log_processing
    assert end_at - start_at < QUEUE_WAIT_SECONDS * 2 + 2

def test_no_cancellation():
    cached_container = os.path.join(os.path.dirname(__file__), "data/logging-example_latest.sif")
    repo_container = "docker://ghcr.io/megosato/logging-example"
    container_uri = cached_container if os.path.exists(cached_container) else repo_container
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri,
            engine_name="singularity",
            kwargs=kwargs,
            file_kwargs={},
            cancel_requested_func=lambda: False,
        )
    }
