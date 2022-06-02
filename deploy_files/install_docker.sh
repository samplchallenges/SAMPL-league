#!/bin/bash -e

# from app/.ebextensions/04_install_docker.config

# 00_update_package_cache
sudo yum update -y

# 10_install_docker
sudo amazon-linux-extras install docker

# 20_start_docker
sudo service docker start

# 30_set_permissions_ec2-user
sudo usermod -a -G docker ec2-user

# 40_set_permissions_webapp
sudo usermod -a -G docker webapp

# start the docker daemon
sudo systemctl start docker 
ps aux | grep docker