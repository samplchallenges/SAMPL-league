#!/bin/bash -e

sudo amazon-linux-extras install nginx1
sudo mv ~/deploy_files/nginx.conf /etc/nginx
sudo nginx -t
sudo systemctl restart nginx
