Admin Tasks
***********

General overview
================

.. figure:: images/overview.png
   :width: 600


When a commit is made to the ``main`` branch (e.g. merging in a pull reqest that targets the ``main`` branch) a build is trigged on `Code Pipeline`_.
Then a notifcation is sent (using `SNS`_) to the ``#sampl-devops`` channel in the Mobley Lab `Slack workspace <https://mobleylab.slack.com>`_.
Once the build completes, the code is sent to the `Elastic Beanstalk`_ instance and the web service restarted.


Amazon Web Services infrastrcture (AWS)
---------------------------------------

All services are hosted in US East (Ohio) us-east-2.


- E-mail: `SES`_
- DNS: `Route 53`_

.. _SNS: https://us-east-2.console.aws.amazon.com/sns/v3/home?region=us-east-2#/dashboard
.. _SES: https://us-east-2.console.aws.amazon.com/sesv2/home?region=us-east-2#/account
.. _Cloud Watch: https://us-east-2.console.aws.amazon.com/cloudwatch/home?region=us-east-2#dashboards:name=Sampl-League
.. _Elastic Beanstalk: https://us-east-2.console.aws.amazon.com/elasticbeanstalk/home?region=us-east-2#/environment/dashboard?applicationName=SAMPL-league&environmentId=e-bkz8t2g9mq
.. _Code Pipeline: https://us-east-2.console.aws.amazon.com/codesuite/codepipeline/pipelines/sampl/view?region=us-east-2
.. _Route 53: https://console.aws.amazon.com/route53/v2/hostedzones#ListRecordSets/Z01835681J808IAHZUIMB


`Python <http://www.python.org/>`_

- `E-mail <https://us-east-2.console.aws.amazon.com/sesv2/home?region=us-east-2#/account>`_
- `Monitering + Alarms <https://us-east-2.console.aws.amazon.com/cloudwatch/home?region=us-east-2#dashboards:name=Sampl-League>`_
- `DNS <https://console.aws.amazon.com/route53/v2/hostedzones#ListRecordSets/Z01835681J808IAHZUIMB>`_
- `Elastic Beanstalk <https://us-east-2.console.aws.amazon.com/elasticbeanstalk/home?region=us-east-2#/environment/dashboard?applicationName=SAMPL-league&environmentId=e-bkz8t2g9mq>`_
- `Code deployment <https://us-east-2.console.aws.amazon.com/codesuite/codepipeline/pipelines/sampl/view?region=us-east-2>`_
- `Server (EC2) <https://us-east-2.console.aws.amazon.com/ec2/v2/home?region=us-east-2#Instances:>`_
- `Object Storage (S3) <https://s3.console.aws.amazon.com/s3/buckets/sampl-league-storage?region=us-east-2&tab=objects>`_
- `Database (RDS) <https://us-east-2.console.aws.amazon.com/rds/home?region=us-east-2#databases:>`_
- `Container Registry (both public and private) <https://us-east-2.console.aws.amazon.com/ecr/repositories?region=us-east-2>`_


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
