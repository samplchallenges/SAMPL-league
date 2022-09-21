#!/bin/bash -e

if [ ! -e /var/app/current/set_staging.sh ]
then
    echo "set_staging.sh must be created before running this command"
fi

sudo -i -u webapp << 'EOF' 
source /var/app/current/.venv/bin/activate

source /var/app/current/set_staging.sh

cd /var/app/current/

python manage.py collectstatic --no-input

python manage.py makemigrations

python manage.py migrate --no-input

EOF
