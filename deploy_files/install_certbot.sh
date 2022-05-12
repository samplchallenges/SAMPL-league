#!/bin/bash -e

# from app/.ebextensions/00_install_certbot.config

# 00_download_epel
wget -r --no-parent -A 'epel-release-*.rpm' http://dl.fedoraproject.org/pub/epel/7/x86_64/Packages/e/

# 10_install_epel_release
rpm -Uvh dl.fedoraproject.org/pub/epel/7/x86_64/Packages/e/epel-release-*.rpm

# 20_enable_epel
yum-config-manager --enable epel*

# 30_install_certbot
yum install -y certbot python2-certbot-nginx
