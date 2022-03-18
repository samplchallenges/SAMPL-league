import os
import subprocess

import pytest

import ever_given.wrapper


ECR_BASE_URL = "103125031445.dkr.ecr.us-east-2.amazonaws.com"
ECR_SAMPLLEAGUE_URL = "103125031445.dkr.ecr.us-east-2.amazonaws.com/sampl-league"


def run_aws_logout(container_engine):
    if container_engine == "singularity":
        os.environ["SINGULARITY_DOCKER_USERNAME"]=""
        os.environ["SINGULARITY_DOCKER_PASSWORD"]=""

def run_aws_login(container_engine):
    if container_engine == "docker":
        subprocess.run(
            [
                f"aws ecr get-login-password --profile sampl_pull | docker login --username AWS --password-stdin {ECR_BASE_URL}"
            ],
            shell=True,
            check=True,
        )
    elif container_engine == "singularity":
        os.environ["SINGULARITY_DOCKER_USERNAME"]="AWS"
        aws_call = subprocess.run(
            ["aws ecr get-login-password --profile sampl_pull"],
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        os.environ["SINGULARITY_DOCKER_PASSWORD"]=aws_call.stdout
    else:
        raise Exception("Container Engine not yet implemented") 


@pytest.mark.parametrize(["container_engine"], [["docker"], ["singularity"]])
def test_awslogin_aws_container(container_engine):
    container_uri = f"{ECR_SAMPLLEAGUE_URL}:calc-molwt_0.0.1"
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri,
            kwargs=kwargs,
            container_type="docker",
            engine_name=container_engine,
            file_kwargs={},
            aws_login_func=run_aws_login
        )
    }
    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(79.04219916)


@pytest.mark.parametrize(["container_engine"], [["docker"], ["singularity"]])
def test_awslogin_ghcr_container(container_engine):
    container_uri = "ghcr.io/megosato/calc-molwt:latest"
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri,
            kwargs=kwargs,
            container_type="docker",
            engine_name=container_engine,
            file_kwargs={},
            aws_login_func=run_aws_logout
        )
    }
    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(79.04219916)




