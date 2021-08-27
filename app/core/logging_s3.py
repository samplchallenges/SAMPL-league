import os.environ
import os.path

import boto3


def upload_log(challenge_id, submission_run_id, evaluation_id, output, error):
    s3 = boto3.client("s3")
    bucket = os.environ["AWS_STORAGE_BUCKET_NAME"]
    key_base = "/logs/{}/{}/{}/".format(challenge_id, submission_run_id, evaluation_id)
    s3.put_object(Body=output, Bucket=bucket, key=os.path.join(key_base, "output.log"))
    s3.put_object(Body=error, Bucket=bucket, key=os.path.join(key_base, "error.log"))
