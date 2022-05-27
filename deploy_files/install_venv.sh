#!/bin/bash -e

sudo -i -u webapp <<'EOF'
/usr/local/bin/pipenv --help
cd /var/app/current
PIPENV_VENV_IN_PROJECT=1 /usr/local/bin/pipenv --rm
PIPENV_VENV_IN_PROJECT=1 /usr/local/bin/pipenv install
source /var/app/current/.venv/bin/activate
pip uninstall -y ever-given
pip install  /var/app/ever_given-0.0.13-py3-none-any.whl
pip install  /var/app/current/dist/sampl_app-0.0.1-py3-none-any.whl
pip install gunicorn
pip install django-debug-toolbar
EOF

source /var/app/current/.venv/bin/activate
pip freeze
