import django
from dask.distributed import Client

django.setup()


def run_submission(submission_id, is_public=True):
    from core.models import (
        Submission,
        Evaluation,
        InputElement,
        InputValue,
        SubmissionRun,
    )
    dask_url = "127.0.0.1:8786"
    client = Client(dask_url)

    submission = Submission.objects.get(pk=submission_id)
    container = submission.container
    challenge = submission.challenge
    if not container.digest:
        container.digest = "DEADBEEF"
        container.save()
    submission_run = SubmissionRun.objects.create(
        submission=submission,
        digest=container.digest,
        is_public=True,
        status=SubmissionRun._Status.PENDING,
    )

    input_data = InputElement.objects.filter(challenge=challenge, is_public=is_public)
    values = [
        InputValue.objects.get(input_element=input_element).value
        for input_element in input_data
    ]
    submissions = []
    for element in input_data:
        future = client.submit(
            run_element, submission, element, submission_run, pure=False
            )
        submissions.append(future)
    # when do we "know" it was a success?
    runs = client.gather(submissions)
    print(runs)
    submission_run._Status.SUCCESS
    submission_run.save()


def run_element(submission, element, submission_run):
    import ever_given.wrapper

    from core.models import (
        FloatValue,
        Prediction,
        Evaluation,
        InputValue,
    )

    container = submission.container
    challenge = submission.challenge
    evaluation = Evaluation.objects.create(
        input_element=element, submission_run=submission_run, exit_status=1
    )
    container_uri = f"{container.registry}/{container.label}:{container.tag}"
    smiles_string = InputValue.objects.get(input_element=element).value
    try:
        command = submission.challenge.execution_options_json["command"]
        command += smiles_string
    except KeyError:  # if no execution options nothing to prepend
        command = smiles_string
    print(command)
    result = ever_given.wrapper.run_submission_container(container_uri, command)
    result = float(result.strip())
    print(result)
    result_obj = FloatValue.objects.create(value=result)
    # maybe we return the prediction so we can score it?
    prediction = Prediction.objects.create(
        challenge=challenge,
        key="molWeight",
        evaluation=evaluation,
        value_object=result_obj,
    )
    evaluation.exit_status = 0
    evaluation.save()
    return 0
