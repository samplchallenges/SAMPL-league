#!/bin/bash -e 

source set_aws_ecr_cred.sh && aws s3 cp s3://aws-container-registry-config/config ~/config && aws s3 cp s3://aws-container-registry-config/credentials ~/credentials
rm set_aws_ecr_cred.sh
sudo mv ~/config ~/credentials /home/webapp/.aws