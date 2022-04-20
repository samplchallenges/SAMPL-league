"""
pip install fabric to run this
"""
from pathlib import Path

from invoke import task
from fabric import Connection


@task
def build(c):
    with c.cd("app"):
        c.run("PIPENV_IGNORE_VIRTUALENVS=1 pipenv run python setup.py bdist_wheel")
    with c.cd("ever_given"):
        c.run("PIPENV_IGNORE_VIRTUALENVS=1 pipenv run python setup.py bdist_wheel")


DEPLOY_DIR = "deploy_files"

PROD_HOST = "sampl.us-east-2.elasticbeanstalk.com"
STAGING_HOST = "ec2-18-191-119-16.us-east-2.compute.amazonaws.com"


def sampl_staging():
    return Connection(
        STAGING_HOST,
        user="ec2-user",
        connect_kwargs={
            "key_filename": "/home/robbason/.ssh/aws-eb",
        },
    )


@task
def deploy(c):
    with sampl_staging() as remote_c:
        remote_c.put("ever_given/dist/ever_given-0.0.13-py3-none-any.whl")
        remote_c.put("app/dist/sampl_app-0.0.1-py3-none-any.whl")
        remote_c.run("ls")


@task
def install_dependencies(c):
    with sampl_staging() as remote_c:
        deploy_files = Path(DEPLOY_DIR)
        remote_c.run(f"mkdir -p {DEPLOY_DIR}")
        for filename in deploy_files.iterdir():
            remote_c.put(filename, remote=str(filename))

        remote_c.run(f"{deploy_files / 'sample_deploy.sh'}")
        remote_c.run("echo 'Hostname:'")
        remote_c.run("hostname")
        remote_c.run("echo 'Files:'")
        remote_c.run("ls")
        remote_c.run("echo 'Can we sudo?'")
        remote_c.run("sudo whoami")
        remote_c.run("sudo amazon-linux-extras install nginx1")
        # remote_c.run("sudo yum install -y nginx")
