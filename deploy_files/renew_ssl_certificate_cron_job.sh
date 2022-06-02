#!/bin/bash -e

# from app/.ebextensions/03_renew_ssl_certificate_cron_job.config

# /tmp/renew_cert_cron
cat > /tmp/renew_cert_cron <<- EOM
0 1,13 * * * certbot renew --no-self-upgrade
EOM

ls -l /tmp/renew_cert_cron

chmod 777 /tmp/renew_cert_cron

# 10_create_cert_crontab
sudo crontab /tmp/renew_cert_cron

# 20_delete_cronjob_file
sudo rm /tmp/renew_cert_cron
