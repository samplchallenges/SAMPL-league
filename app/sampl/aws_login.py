import os.environ
import subprocess


def run_aws_login(ecr_base_url, container_engine):
    if container_engine == "docker":
        subprocess.run(
            [
                f"aws ecr get-login-password --region us-east-2 --profile sampl_pull | docker login --username AWS --password-stdin {ecr_base_url}"
            ],
            shell=True,
            check=True,
        )
    elif container_engine == "singularity":
        os.environ["SINGULARITY_DOCKER_USERNAME"] = "AWS"
        aws_call = subprocess.run(
            [
                "aws",
                "ecr",
                "get-login-password",
                "--region",
                "us-east-2",
                "--profile",
                "sampl_pull",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        os.environ["SINGULARITY_DOCKER_PASSWORD"] = aws_call.stdout
    else:
        raise Exception("Container Engine not yet implemented")


def run_aws_logout(container_engine):
    if container_engine == "singularity":
        os.environ["SINGULARITY_DOCKER_USERNAME"] = ""
        os.environ["SINGULARITY_DOCKER_PASSWORD"] = ""
