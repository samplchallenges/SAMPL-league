#!/bin/bash -e 

# 01_create_file
sudo touch /var/log/app.log

# 02_change_permissions
sudo chmod g+s /var/log/app.log

# 03_change_owner
sudo chown webapp:webapp /var/log/app.log
