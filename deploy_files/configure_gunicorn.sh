#!/bin/bash -e

sudo mkdir /var/pids
sudo chmod 777 /var/pids

sudo mkdir /opt/env
sudo mv ~/deploy_files/env /opt/env/env
cat /opt/env/env
sudo mv ~/deploy_files/gunicorn.socket /etc/systemd/system
sudo mv ~/deploy_files/gunicorn.service /etc/systemd/system
sudo mv ~/deploy_files/jupyter.service /etc/systemd/system
sudo mv ~/deploy_files/scheduler.service /etc/systemd/system
sudo mv ~/deploy_files/worker.service /etc/systemd/system

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
