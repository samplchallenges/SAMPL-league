#!/bin/bash -e

# from app/.ebextensions/00_install_certbot.config

# 00_download_epel
sudo wget -r --no-parent -A 'epel-release-*.rpm' http://dl.fedoraproject.org/pub/epel/7/x86_64/Packages/e/

# 10_install_epel_release
sudo rpm -Uvh dl.fedoraproject.org/pub/epel/7/x86_64/Packages/e/epel-release-*.rpm

# 20_enable_epel
sudo yum-config-manager --enable epel*

# 30_install_certbot
sudo yum install -y certbot python2-certbot-nginx
