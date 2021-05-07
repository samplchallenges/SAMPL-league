import django
import dask.bag
from dask.distributed import Client, worker_client


import ever_given.wrapper

from core.models import (
    AnswerKey,
    Evaluation,
    FloatValue,
    InputElement,
    InputValue,
    Prediction,
    Submission,
    SubmissionRun,
)

# should this happen on the module level or in the function/method?
# django.setup()


def score_submission(submission_id, *run_ids):

    submission = Submission.objects.get(pk=submission_id)
    challenge = submission.challenge
    smiles_type = challenge.inputtype_set.get(key="SMILES")
    for run_id in run_ids:

        submission_run = submission.submissionrun_set.get(pk=run_id)
        evaluations = submission_run.evaluation_set.all()
        for evaluation in evaluations:
            prediction = evaluation.prediction_set.get(key="molWeight")
            input_element = evaluation.input_element
            input_value = input_element.inputvalue_set.get(input_type=smiles_type)

            answer_key = challenge.answerkey_set.get(
                input_element=input_element, key="molWeight"
            )
            print(
                "Prediction:",
                prediction.value,
                input_value.value,
                "answer:",
                answer_key.value,
            )


def get_status(key):
    dask_url = "127.0.0.1:8786"
    client = Client(dask_url)
    return client.get_metadata(key)


def run_submission(submission_id, is_public=True):

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
        is_public=is_public,
        status=SubmissionRun._Status.PENDING,
    )
    input_elements = challenge.inputelement_set.filter(
        is_public=submission_run.is_public
    )
    args = [
        (submission.id, input_element.id, submission_run.id)
        for input_element in input_elements
    ]
    elements_bag = dask.bag.from_sequence(args, partition_size=1)
    dask_graph = elements_bag.map(run_element)
    dask_graph.visualize(filename="bag_graph.svg")
    print("submitting task")
    future = client.compute(dask_graph)
    # client.submit(
    # run_all_elements, submission.id, submission_run.id, pure=False
    # )

    submission_run.digest = future.key
    submission_run.save()
    return submission_run


def run_all_elements(submission_id, submission_run_id):

    submission = Submission.objects.get(pk=submission_id)
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)
    input_elements = challenge.inputelement_set.filter(
        is_public=submission_run.is_public
    )
    args = [
        (submission.id, input_element.id, submission_run.id)
        for input_element in input_elements
    ]
    print("about to submit!")
    with worker_client() as client:
        futures = client.map(run_element, args, pure=False)

        print("futures", futures)
        print("About to gather")
        # when do we "know" it was a success?
        runs = client.gather(futures)
    submission_run.status = SubmissionRun._Status.SUCCESS
    submission_run.save()
    print(runs)


def run_element(submission_id, element_id, submission_run_id):

    submission = Submission.objects.get(pk=submission_id)
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)
    element = challenge.inputelement_set.get(
        pk=element_id, is_public=submission_run.is_public
    )

    container = submission.container

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
