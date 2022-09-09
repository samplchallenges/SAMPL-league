import os.path
from pathlib import Path
import tempfile
import subprocess
import shutil

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
    tempdir = Path(tempfile.mkdtemp())
    container_sif_path = tempdir / "calc-molwt_latest-unique.sif"
    command = ["singularity", "pull", container_sif_path, f"docker://{container_uri}"]
    subprocess.run(command, check=True)
    assert os.path.exists(container_sif_path)
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_sif_path,
            kwargs=kwargs,
            container_type="singularity_local",
            engine_name="singularity",
            file_kwargs={}
        )
    }
    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(79.04219916)
    shutil.rmtree(tempdir)


def test_pull_container_singularity_engine():
    container_uri = "ghcr.io/megosato/calc-molwt:latest" 
    tempdir = Path(tempfile.mkdtemp())
    container_sif_path = tempdir / "calc-molwt_latest-unique.sif"
    ever_given.wrapper.pull_container(container_uri, "docker", "singularity", container_sif_path)

    assert os.path.exists(container_sif_path)

    shutil.rmtree(tempdir)


def test_pull_container_docker_engine():
    container_uri = "ghcr.io/megosato/calc-molwt:latest" 
    ever_given.wrapper.pull_container(container_uri, "docker", "docker", None)

    command = ['docker', 'pull', container_uri]
    ended_proc = subprocess.run(command, capture_output=True, check=True)

    output = ended_proc.stdout.decode('utf-8')
    assert f'Status: Image is up to date for {container_uri}' in output


