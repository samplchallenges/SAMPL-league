from dask.distributed import Client
import sys

def fire_off_tasks(submission_id, elements):
    dask_url = "127.0.0.1:8786"
    client = Client(dask_url)
    #client.upload_file("code.zip")
    #client.upload_file("referee/tasks.py")
    client.upload_file("tasks.py")
    submissions = []
    for element in elements:
        future = client.submit(
            run_submission_element, submission_id, element, pure=False
        )
        submissions.append(future)
    print("gathering")
    print(sys.path)
    results = client.gather(submissions)
    print(results)
    print("done")


def run_submission_element(submission_id, element):
    print("fired")
    from core.models import Submission
    #import ever_given

    print("running")
    submission = Submission.objects.get(pk=submission_id)
    container = submission.container
    container_uri = f"{container.registry}/{container.label}:{container.tag}"
    command = submission.challenge.execution_options_json["command"]
    print(command)
    print(container_uri)
    result = ever_given.run_submission_container(container_uri, command)
    return result
