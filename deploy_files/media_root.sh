#!/bin/bash -e

# from: app/.ebextensions/06_media_root.config

# 01_create_media_dir
sudo mkdir -p /opt/app/sampl/media

# 02_set_media_dir_ownership:
sudo chown webapp:webapp /opt/app/sampl/media
