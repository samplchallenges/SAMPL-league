#!/bin/bash -e

sudo amazon-linux-extras install python3.8
sudo yum install python3-pip -y
#sudo yum reinstall python3-pip -y
sudo -H pip3 install -U pipenv
pipenv
