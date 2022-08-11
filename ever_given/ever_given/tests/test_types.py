import os.path
import tempfile
import subprocess

import pytest

import ever_given.wrapper


def test_docker_container(container_engine):
    container_uri = "ghcr.io/megosato/calc-molwt:latest"
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri, 
            kwargs=kwargs,
            container_type="docker",
            engine_name=container_engine,
            file_kwargs={}
        )
    }
    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(79.04219916)


def test_singularity_sif_container_docker_engine():
    container_sif = "calc-molwt_latest.sif"
    kwargs = {"smiles": "c1cccnc1"}
    with pytest.raises(ValueError):
        results = {
            key: value
            for key, value in ever_given.wrapper.run(
                container_sif,
                kwargs=kwargs,
                container_type="singularity_local",
                engine_name="docker",
                file_kwargs={}
            )
        }


def test_singularity_sif_container_singularity_engine():
    container_uri = "ghcr.io/megosato/calc-molwt:latest" 
    container_sif = "calc-molwt_latest-unique.sif"
    command = ["singularity", "pull", container_sif, f"docker://{container_uri}"]
    subprocess.run(command, check=True)
    assert os.path.exists(container_sif)
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_sif,
            kwargs=kwargs,
            container_type="singularity_local",
            engine_name="singularity",
            file_kwargs={}
        )
    }
    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(79.04219916)
