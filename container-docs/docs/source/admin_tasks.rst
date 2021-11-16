Admin Tasks
***********

# General overview

.. figure:: images/overview.png
   :width: 500

# AWS infra, when needed region is US East (Ohio) us-east-2

E-mail https://us-east-2.console.aws.amazon.com/sesv2/home?region=us-east-2#/account
Monitering + Alarms https://us-east-2.console.aws.amazon.com/cloudwatch/home?region=us-east-2#dashboards:name=Sampl-League
DNS https://console.aws.amazon.com/route53/v2/hostedzones#ListRecordSets/Z01835681J808IAHZUIMB
Elastic Beanstalk https://us-east-2.console.aws.amazon.com/elasticbeanstalk/home?region=us-east-2#/environment/dashboard?applicationName=SAMPL-league&environmentId=e-bkz8t2g9mq
Code deployment https://us-east-2.console.aws.amazon.com/codesuite/codepipeline/pipelines/sampl/view?region=us-east-2
Server (EC2) https://us-east-2.console.aws.amazon.com/ec2/v2/home?region=us-east-2#Instances:
Object Storage (S3) https://s3.console.aws.amazon.com/s3/buckets/sampl-league-storage?region=us-east-2&tab=objects
Database (RDS) https://us-east-2.console.aws.amazon.com/rds/home?region=us-east-2#databases:
Container Registry (both public and private) https://us-east-2.console.aws.amazon.com/ecr/repositories?region=us-east-2


# Update ever_given

cd ever_given
# make changes needed for PR
# bump version
then 
python -m pip install --upgrade build
python -m build
python3 -m twine upload  dist/*

You will need an account on https://pypi.org/
The new package should show up here https://pypi.org/project/ever-given/#history

# Spin up worker


# Expand local storage
Adapted from here
https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/recognize-expanded-volume-linux.html


df -hT
lsblk

sudo growpart /dev/nvme0n1 1

df -hT
lsblk

sudo xfs_growfs -d /
