#!/bin/bash -e

sudo mv /home/ec2-user/deploy_files/set_staging.sh /var/app/current/

sudo -i -u webapp << 'EOF' 
source /var/app/current/.venv/bin/activate

source /var/app/current/set_staging.sh

cd /var/app/current/

python manage.py createsuperuser --no-input

python manage.py sample_data --delete
EOF
