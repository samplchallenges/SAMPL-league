"""
pip install fabric to run this
"""
from pathlib import Path

from invoke import task
from fabric import Connection
from patchwork.files import exists
from patchwork.transfers import rsync

@task
def build(c):
    with c.cd("app"):
        c.run("PIPENV_IGNORE_VIRTUALENVS=1 pipenv run python setup.py bdist_wheel")
        c.run("tar -cf ../current.tar manage.py Pipfile Pipfile.lock Procfile daskworkerinit.py .platform dist")
    with c.cd("ever_given"):
        c.run("PIPENV_IGNORE_VIRTUALENVS=1 pipenv run python setup.py bdist_wheel")

DEPLOY_DIR = "deploy_files"

PROD_HOST = "sampl.us-east-2.elasticbeanstalk.com"
STAGING_HOST = "ec2-3-134-189-93.us-east-2.compute.amazonaws.com"

KEY_FILE =  "/Users/megosato/.ssh/aws-eb"


def sampl_staging(user="ec2-user"):
    return Connection(
        STAGING_HOST,
        user=user,
        connect_kwargs={
        "key_filename": KEY_FILE,
        },
    )


@task
def deploy(c):
    with sampl_staging() as remote_c:
        remote_c.put("ever_given/dist/ever_given-0.0.13-py3-none-any.whl")
        if not exists(remote_c, "/var/app"):
            remote_c.sudo("mkdir /var/app/")
        if not exists(remote_c, "/var/app/current"):
            remote_c.sudo("mkdir /var/app/current")
        remote_c.put("current.tar")
        remote_c.sudo("mv current.tar /var/app/current")
        remote_c.sudo("mv ever_given-0.0.13-py3-none-any.whl /var/app")
        with remote_c.cd("/var/app/current"):
            remote_c.run("sudo tar -xvf current.tar")
            remote_c.run("sudo rm current.tar")
            remote_c.run("ls -l .")


# def _install_venv(remote_c):
#     with remote_c.cd("/var/app/current"):
#         remote_c.run("pipenv install")
#         remote_c.run("pipenv run pip uninstall ever-given")
#         remote_c.run("pipenv run pip install ../ever_given-0.0.13-py3-none-any.whl")
#         remote_c.run("pipenv run pip install dist/sampl_app-0.0.1-py3-none-any.whl")
#         remote_c.run("pipenv run pip freeze")

def _download_aws_ecr_cred(c):
    with sampl_staging() as remote_c:
        if not exists(remote_c, "/home/webapp/.aws", sudo=True):
            remote_c.sudo("mkdir /home/webapp/.aws")
        remote_c.put("app/set_aws_ecr_cred.sh")
        remote_c.run("source set_aws_ecr_cred.sh && aws s3 cp s3://aws-container-registry-config/config ~/config && aws s3 cp s3://aws-container-registry-config/credentials ~/credentials")
        remote_c.run("rm set_aws_ecr_cred.sh")
        remote_c.run("sudo mv ~/config ~/credentials /home/webapp/.aws")

def _set_env_variables(c):
    with sampl_staging() as remote_c:
        remote_c.put("app/set_staging.sh")
        remote_c.run("echo source /home/ec2-user/set_staging.sh >> ~/.bashrc")
        remote_c.run("source ~/.bashrc")

def _upload_dependency_install_file(remote_c, filename):
    if not exists(remote_c, "deploy_files"):
        remote_c.run("mkdir deploy_files")
    filepath = Path("deploy_files") / filename
    remote_c.put(filepath, remote=str(filepath))

def _run_install_file(remote_c, sh_file, sudo=False):
    deploy_files = Path("deploy_files")
    cmd = f"sudo bash {deploy_files / sh_file}" if sudo else f"bash {deploy_files / sh_file}" 
    remote_c.run(cmd)

def _install_dependency(install_file, sudo=False):
    with sampl_staging() as remote_c:
        _upload_dependency_install_file(remote_c, install_file)
        _run_install_file(remote_c, install_file, sudo)

@task
def install_singularity(c):
    install_file = "install_singularity.sh"
    _install_dependency(install_file)

@task 
def install_python38(c):
    install_file = "install_python38.sh"
    _install_dependency(install_file)

@task
def install_docker(c):
    install_file = "install_docker.sh"
    _install_dependency(install_file)

@task
def install_certbot(c):
    install_file = "logging.sh"
    _install_dependency(install_file)

@task
def open_https_port(c):
    install_file = "open_https_port.sh"
    _install_dependency(install_file)

@task
def grant_executable_rights(c):
    install_file = "grant_executable_rights.sh"
    _install_dependency(install_file)

@task
def renew_ssl_certificate_cron_job(c):
    install_file = "renew_ssl_certificate_cron_job.sh"
    _install_dependency(install_file)

@task
def setup_logging(c):
    install_file = "logging.sh"
    _install_dependency(install_file)

@task
def setup_media_root(c):
    install_file = "media_root.sh"
    _install_dependency(install_file)

def install_venv(c):
    install_file = "install_venv.sh"
    _install_dependency(install_file)

@task
def install_dependencies(c):
    with sampl_staging() as remote_c:
        _set_env_variables(c)
        setup_logging(c)
        install_venv(c)
    #_download_aws_ecr_cred(c)
