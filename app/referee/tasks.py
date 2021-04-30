from dask.distributed import Client


def fire_off_tasks(submission_id, elements):
    dask_url = "127.0.0.1:8786"
    client = Client(dask_url)
    submissions = []
    for element in elements:
        future = client.submit(
            run_submission_element, submission_id, element, pure=False
        )
        submissions.append(future)
    print("gathering")
    results = client.gather(submissions)
    print(results)
    print("done")


def run_submission_element(submission_id, element):
    print("fired")
    import ever_given.wrapper

    from core.models import (
        Evaluation,
        FloatValue,
        InputElement,
        Prediction,
        Submission,
        SubmissionRun,
    )

    print("running")

    submission = Submission.objects.get(pk=submission_id)
    container = submission.container
    challenge = submission.challenge
    # This should be passed in as the element
    input_element = InputElement.objects.get(
        name="benzene",
        challenge=challenge
    )

    if not container.digest:
        container.digest = "DEADBEEF"
        container.save()
    submission_run = SubmissionRun.objects.create(
        submission=submission,
        digest=container.digest,
        is_public=True,
        status=SubmissionRun._Status.PENDING,
    )

    evaluation = Evaluation.objects.create(
        input_element=input_element, submission_run=submission_run, exit_status=0
    )

    container_uri = f"{container.registry}/{container.label}:{container.tag}"
    try:
        command = submission.challenge.execution_options_json["command"]
        command += element
    except KeyError:  # if no execution options nothing to prepend
        command = element
    print(command)
    print(container_uri)
    result = ever_given.wrapper.run_submission_container(container_uri, command)
    result = float(result.strip())
    result_obj = FloatValue.objects.create(value=result)
    #prediction = Prediction.objects.create(
    #    challenge=challenge,
    #    key="molWeight",
    #    content_type=FloatValue,
    #    evaluation=evaluation,
    #)

    submission_run._Status.SUCCESS
    submission_run.save()
    return (element, result)
