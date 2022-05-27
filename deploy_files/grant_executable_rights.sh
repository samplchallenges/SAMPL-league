#!/bin/bash -e

# from app/.ebextensions/02_grant_executable_rights.config
sudo chmod +x /var/app/current/.platform/hooks/postdeploy/00_get_cert.sh
