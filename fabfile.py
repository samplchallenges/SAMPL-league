"""
pip install fabric to run this
"""
import os
import os.path
from pathlib import Path
import tempfile

from invoke import task
from fabric import Connection
from patchwork.files import exists
from patchwork.transfers import rsync

@task
def build(c):
    with c.cd("app"):
        tar_cmd = """tar -cf ../current.tar \
                        Pipfile Pipfile.lock Procfile \
                        core notebooks referee sampl dist .platform \
                        daskworkerinit.py manage.py setup.py"""
        c.run("PIPENV_IGNORE_VIRTUALENVS=1 pipenv run python setup.py bdist_wheel")
        c.run(tar_cmd)
    with c.cd("ever_given"):
        c.run("PIPENV_IGNORE_VIRTUALENVS=1 pipenv run python setup.py bdist_wheel")

DEPLOY_DIR = "deploy_files"

PROD_HOST = "sampl.us-east-2.elasticbeanstalk.com"
STAGING_HOST = "app-staging.samplchallenges.org"

HOME = os.environ.get("HOME")
default_key_file = os.path.join(HOME, ".ssh/aws-eb")
KEY_FILE = os.environ.get("SAMPL_SSH_KEY_FILE", default_key_file)

def sampl_staging(user="ec2-user"):
    return Connection(
        STAGING_HOST,
        user=user,
        connect_kwargs={
        "key_filename": KEY_FILE,
        },
    )


# TODO: pick up ever-given version automatically
@task
def deploy(c):
    with sampl_staging() as remote_c:
        remote_c.put("ever_given/dist/ever_given-0.0.14-py3-none-any.whl")
        if not exists(remote_c, "/var/app"):
            remote_c.sudo("mkdir /var/app/")
        if not exists(remote_c, "/var/app/current"):
            remote_c.sudo("mkdir /var/app/current")
        remote_c.put("current.tar")
        remote_c.sudo("mv current.tar /var/app/current")
        remote_c.sudo("mv ever_given-0.0.14-py3-none-any.whl /var/app")
        with remote_c.cd("/var/app/current"):
            remote_c.run("sudo tar -xvf current.tar")
            remote_c.run("sudo rm current.tar")
            remote_c.run("ls -l .")
        with remote_c.cd("/var/app/"):
            remote_c.run("sudo chown -R webapp:webapp current")

def _upload_dependency_install_file(remote_c, filename):
    if not exists(remote_c, "deploy_files"):
        remote_c.run("mkdir deploy_files")
    filepath = Path("deploy_files") / filename
    remote_c.put(filepath, remote=str(filepath))

def _run_install_file(remote_c, sh_file, *args, sudo=False):
    deploy_files = Path("deploy_files")
    sh_file = f"{deploy_files / sh_file}"
    cmd = f"sudo {sh_file}" if sudo else sh_file
     # TODO: join args in a more shell-friendly way if args may need it
    full_cmd = cmd + " " + " ".join(args)
    try: 
        remote_run = remote_c.run(full_cmd)
    except:
        raise Exception(f"DeployError: {sh_file}, {full_cmd}")

def _install_dependency(install_file, *args, sudo=False):
    with sampl_staging() as remote_c:
        _upload_dependency_install_file(remote_c, install_file)
        import pdb;pdb.set_trace()
        _run_install_file(remote_c, install_file, *args, sudo=sudo)

@task
def download_aws_ecr_credentials(c):
    with sampl_staging() as remote_c:
        if not exists(remote_c, "/home/webapp/.aws", sudo=True):
            remote_c.sudo("mkdir /home/webapp/.aws")
        remote_c.put("app/set_aws_ecr_cred.sh")
    install_file = "download_aws_ecr_credentials.sh"
    _install_dependency(install_file)

@task
def create_webapp_user(c):
    with sampl_staging() as remote_c:
        install_file = "create_webapp_user.sh"
        _install_dependency(install_file)

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
    install_file = "install_certbot.sh"
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

@task
def install_venv(c):
    install_file = "install_venv.sh"
    _install_dependency(install_file)

@task
def setup_djangoapp(c):
    install_file = "setup_djangoapp.sh"
    _install_dependency(install_file)


@task
def django_initial_data(c):
    with sampl_staging() as remote_c:
        _upload_dependency_install_file(remote_c, "set_staging.sh")
    install_file = "django_initial_data.sh"
    _install_dependency(install_file)


@task
def configure_nginx(c):
    with sampl_staging() as remote_c:
        _upload_dependency_install_file(remote_c, "nginx.conf")
    install_file = "configure_nginx.sh"
    _install_dependency(install_file)

@task
def configure_gunicorn(c):
    with sampl_staging() as remote_c:
        _upload_dependency_install_file(remote_c, "env")
        _upload_dependency_install_file(remote_c, "gunicorn.socket")
        _upload_dependency_install_file(remote_c, "gunicorn.service")
        _upload_dependency_install_file(remote_c, "jupyter.service")
        _upload_dependency_install_file(remote_c, "scheduler.service")
        _upload_dependency_install_file(remote_c, "worker.service")
    install_file = "configure_gunicorn.sh"
    _install_dependency(install_file)


@task
def get_cert(c):
    install_file = "get_cert.sh"
    _install_dependency(install_file)

@task
def restart_all(c):
    # with sampl_staging() as remote_c:
    #     _upload_dependency_install_file(remote_c, "env")
    install_file = "restart_processes.sh"
    _install_dependency(install_file)

@task
def install_dependencies(c):
    install_certbot(c)                  #00_install_certbot.config
    open_https_port(c)                  #01_open_https_port.config
    grant_executable_rights(c)          #02_grant_executable_rights.config
    renew_ssl_certificate_cron_job(c)   #03_renew_ssl_certificate_cron_job.config
    install_docker(c)                   #04_install_docker.config
    install_singularity(c)              #05_install_singularity.config
    setup_logging(c)                    #06_logging.config
    setup_media_root(c)                 #07_media_root.config

    download_aws_ecr_credentials(c)     

    install_python38(c)
    install_venv(c)
    setup_djangoapp(c)

    configure_nginx(c)
    configure_gunicorn(c)
    get_cert(c)


@task
def deploy_full_pipeline(c):
    create_webapp_user(c)
    build(c)
    deploy(c)
    install_dependencies(c)
    django_initial_data(c)


@task
def redeploy_pipeline_venv(c):
    build(c)
    deploy(c)
    install_venv(c)
    setup_djangoapp(c)
    restart_all(c)


@task(pre=[build, deploy, setup_djangoapp, restart_all])
def redeploy_samplapp(c):
    print("Redeployed")


@task
def deploy_env_var(c):
    restart_all(c)


@task
def publish_challenge(c, config_yml):
    """Uploads the directory containing the config_yml"""
    config_dir = Path(config_yml).parent
    config_basename = Path(config_yml).name
    with sampl_staging() as remote_c:
        _upload_dependency_install_file(remote_c, "load_yaml.sh")
        remote_name = "/tmp/challenge_config.tar.gz"
        with tempfile.NamedTemporaryFile(suffix=".tar.gz") as temp_tarfile:
            with c.cd(config_dir):
                c.run(f"tar -czf {temp_tarfile.name} .")
            remote_c.put(temp_tarfile.name, remote_name)
        # TODO: could clash if more than one using same path at the same time
        remote_dir = "/tmp/challenge_config"
        #remote_c.run("rm -rf {remote_dir}")
        remote_c.run(f"mkdir -p {remote_dir}")
        with remote_c.cd(remote_dir):
            remote_c.run(f"tar xzf {remote_name}")

        cmd = f"deploy_files/load_yaml.sh {remote_dir}/{config_basename}"
        remote_c.run(cmd)