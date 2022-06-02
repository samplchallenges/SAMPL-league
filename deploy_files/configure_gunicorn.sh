#!/bin/bash -e

DEPLOY_FILES_DIR="/home/ec2-user/deploy_files"

if ls /var/pids; then
   :
else
   sudo mkdir /var/pids
   sudo chmod 777 /var/pids
fi

if ls /opt/env; then
   :
else
   sudo mkdir /opt/env
fi
sudo mv $DEPLOY_FILES_DIR/env /opt/env/env
cat /opt/env/env
sudo mv $DEPLOY_FILES_DIR/gunicorn.socket /etc/systemd/system
sudo mv $DEPLOY_FILES_DIR/gunicorn.service /etc/systemd/system
sudo mv $DEPLOY_FILES_DIR/jupyter.service /etc/systemd/system
sudo mv $DEPLOY_FILES_DIR/scheduler.service /etc/systemd/system
sudo mv $DEPLOY_FILES_DIR/worker.service /etc/systemd/system

sudo systemctl daemon-reload

sudo systemctl enable --now gunicorn.service
ps aux | grep gunicorn
sleep 3
sudo systemctl enable --now jupyter.service
ps aux | grep jupyter
sleep 3
sudo systemctl enable --now scheduler.service
ps aux | grep dask-scheduler
sleep 3
sudo systemctl enable --now worker.service
ps aux | grep dask-worker
