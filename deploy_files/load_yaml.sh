#!/bin/bash -e

if [ ! -e /var/app/current/set_staging.sh ]
then
    echo "set_staging.sh must be created before running this command"
fi

sudo -i -u webapp << EOF
source /var/app/current/.venv/bin/activate

source /var/app/current/set_staging.sh

cd /var/app/current/
export DJANGO_SETTINGS_MODULE=sampl.settings_prod_console_log  # Otherwise no logging to command line
python manage.py load_yaml --delete $1

EOF