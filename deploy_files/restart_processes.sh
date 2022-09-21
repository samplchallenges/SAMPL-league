#!/bin/bash -e

# DEPLOY_FILES_DIR="/home/ec2-user/deploy_files"

# sudo mv $DEPLOY_FILES_DIR/env /opt/env/env

#sudo pkill -HUP gunicorn

sudo systemctl daemon-reload
sudo systemctl restart gunicorn.service
sudo systemctl restart worker.service
sudo systemctl restart scheduler.service
sudo systemctl restart jupyter.service
sudo systemctl restart docker
