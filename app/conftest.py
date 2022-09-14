import os
from collections import namedtuple
from datetime import datetime, timezone

import dask.distributed as dd
import pytest
from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
from django.core.files.uploadedfile import SimpleUploadedFile

from core import models

# Collection tuple to simplify test setup
TestConfig = namedtuple(
    "TestConfig", ["challenge", "input_type", "output_type", "submission_run"]
)


@pytest.fixture(scope="session")
def dask_client():
    # We need at least 4 workers
    with dd.LocalCluster(n_workers=4, preload=("daskworkerinit_tst.py",)) as cluster:
        yield dd.Client(cluster)


@pytest.fixture
def users(db):
    User = get_user_model()
    user, _ = User.objects.get_or_create(username="hello")
    other_user, _ = User.objects.get_or_create(username="other")
    return user, other_user


@pytest.fixture
def user(users):
    return users[0]


@pytest.fixture
def other_user(users):
    return users[1]


@pytest.fixture
def challenge(challenge_factory, db):
    return challenge_factory("SAMPL1")


@pytest.fixture
def challenge_factory(db):
    def maker(name):
        empty_url = "http://github.com"
        start_at = datetime(2020, 1, 1, hour=1, tzinfo=timezone.utc)
        end_at = datetime(2020, 9, 1, hour=1, tzinfo=timezone.utc)
        challenge = models.Challenge(
            name=name,
            start_at=start_at,
            end_at=end_at,
            repo_url=empty_url,
        )
        challenge.save()
        return challenge

    return maker


@pytest.fixture
def config_factory(challenge_factory, container_factory, db):
    """
    Create a challenge and related objects for testing
    """

    def maker(
        challenge_name,
        label,
        score_label,
        input_key,
        input_model,
        batch_method,
        output_key,
        output_model,
    ):
        challenge = challenge_factory(challenge_name)
        scoring_container = container_factory(challenge, score_label, tag="latest")
        models.ScoreMaker.objects.create(
            challenge=challenge, container=scoring_container
        )
        for key, level in (
            ("diff", models.ScoreType.Level.EVALUATION),
            ("rmse", models.ScoreType.Level.SUBMISSION_RUN),
        ):
            models.ScoreType.objects.create(challenge=challenge, key=key, level=level)

        input_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=True,
            content_type=ContentType.objects.get_for_model(input_model),
            key=input_key,
            description="N/A",
            batch_method=batch_method,
        )
        output_type = models.ValueType.objects.create(
            challenge=challenge,
            is_input_flag=False,
            content_type=ContentType.objects.get_for_model(output_model),
            key=output_key,
            description="N/A",
        )
        container = container_factory(challenge, label, tag="latest")
        submission = models.Submission.objects.create(
            name="Draft Submission",
            user=container.user,
            container=container,
            challenge=challenge,
        )
        submission_run = models.SubmissionRun.objects.create(
            submission=submission,
            digest="cafef00d",
            is_public=True,
            status=models.Status.PENDING,
        )
        return TestConfig(challenge, input_type, output_type, submission_run)

    return maker


@pytest.fixture
def container_factory(user, db):
    def container_maker(challenge, label, tag):
        return models.Container.objects.create(
            name="Container1",
            user=user,
            challenge=challenge,
            container_type=models.ContainerType.DOCKER,
            registry="ghcr.io",
            label=label,
            tag=tag,
        )

    return container_maker


@pytest.fixture
def container(container_factory, challenge):
    return container_factory(
        challenge,
        label="robbason/calc-molwt",
        tag="latest",
    )


@pytest.fixture
def container_arg_factory(db):
    def container_arg_maker(
        container, key, string_value=None, file_name=None, file_body=None
    ):
        if file_name:
            file_value = SimpleUploadedFile(file_name, file_body.encode())
        else:
            file_value = None
        return models.ContainerArg.objects.create(
            container=container,
            key=key,
            file_value=file_value,
            string_value=string_value,
        )

    return container_arg_maker


@pytest.fixture
def scoring_container(challenge, user, db):
    return models.Container.objects.create(
        name="subtraction container",
        user=user,
        container_type=models.ContainerType.DOCKER,
        challenge=challenge,
        registry="ghcr.io",
        label="megosato/score-coords",
        tag="latest",
    )


@pytest.fixture
def score_maker(challenge, scoring_container, db):
    return models.ScoreMaker.objects.create(
        challenge=challenge, container=scoring_container
    )


@pytest.fixture
def draft_submission(submission_factory, container, db):
    return submission_factory(container)


@pytest.fixture
def molfile_type(challenge, db):
    return models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=True,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="molfile",
        description="SD file",
    )


@pytest.fixture
def testing_data_path():
    return os.path.join(os.path.dirname(__file__), "core", "tests", "data")


@pytest.fixture
def elem_factory(testing_data_path, db):
    def elem_maker(
        challenge,
        file_type,
        name,
        file_name,
        parent=None,
        is_parent=False,
        is_public=True,
    ):
        elem = models.InputElement.objects.create(
            name=name,
            challenge=challenge,
            is_public=is_public,
            is_parent=is_parent,
            parent=parent,
        )
        elem.full_clean()
        file_path = os.path.join(testing_data_path, file_name)
        file_value = models.FileValue.from_string(
            file_path, challenge=challenge, input_element=elem
        )
        file_value.save()
        models.InputValue.objects.create(
            input_element=elem, value_type=file_type, value_object=file_value
        )
        return elem

    return elem_maker


@pytest.fixture
def submission_factory(db):
    def submission_maker(container):
        submission = models.Submission(
            name="Draft Submission",
            user=container.user,
            container=container,
            challenge=container.challenge,
        )
        submission.full_clean()
        submission.save()
        return submission

    return submission_maker


@pytest.fixture
def submission_run_factory(db):
    def srun_maker(submission):
        return models.SubmissionRun.objects.create(
            submission=submission,
            digest="cafef00d",
            is_public=True,
            status=models.Status.PENDING,
        )

    return srun_maker


@pytest.fixture
def float_answer_key_factory(db):
    def fak_maker(challenge, elem, value_type, value):
        float_value = models.FloatValue.from_string(
            "72.0", challenge=challenge, input_element=elem
        )
        float_value.save()
        answer_key = models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=value_type,
            value_object=float_value,
        )
        return answer_key

    return fak_maker


@pytest.fixture
def file_answer_key_factory(testing_data_path, db):
    def fak_maker(challenge, elem, value_type, file_name):
        file_path = os.path.join(testing_data_path, file_name)
        file_value = models.FileValue.from_string(
            file_path, challenge=challenge, input_element=elem
        )
        answer_key = models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=value_type,
            value_object=file_value,
        )
        return answer_key

    return fak_maker


@pytest.fixture
def molfile_molw_config(config_factory):
    return config_factory(
        "molfile_molw",
        "robbason/calc-molwt",
        "robbason/score-coords",
        "molfile",
        models.FileValue,
        "sdf",
        "molWeight",
        models.FloatValue,
    )


@pytest.fixture
def smiles_molw_config(config_factory):
    return config_factory(
        "smiles_molw",
        "robbason/calc-molwt",
        "robbason/score-coords",
        "smiles",
        models.TextValue,
        "csv",
        "molWeight",
        models.FloatValue,
    )


@pytest.fixture
def smiles_docking_config_and_func(config_factory, elem_factory):
    config = config_factory(
        "smiles_docking",
        "robbason/calc-molwt",
        "megosato/score-coords",
        "smiles",
        models.TextValue,
        "csv",
        "molWeight",
        models.FloatValue,
    )
    molw_type = config.output_type

    protein_type = models.ValueType.objects.create(
        challenge=config.challenge,
        is_input_flag=True,
        on_parent_flag=True,
        content_type=ContentType.objects.get_for_model(models.FileValue),
        key="protein_pdb",
        description="Protein structure (shared among input elements)",
    )
    parent = elem_factory(
        config.challenge,
        protein_type,
        "5QCR",
        "5qcr.pdb",
        is_parent=True,
        is_public=True,
    )

    def add_element(name, smiles, mol_weight):
        element = models.InputElement.objects.create(
            challenge=config.challenge,
            parent=parent,
            is_parent=False,
            name=name,
            is_public=True,
        )
        smiles_value = models.TextValue.from_string(
            smiles, challenge=config.challenge, input_element=element
        )
        smiles_value.save()
        models.InputValue.objects.create(
            input_element=element,
            value_type=config.input_type,
            value_object=smiles_value,
        )
        float_value = models.FloatValue.from_string(
            str(mol_weight), challenge=config.challenge, input_element=element
        )
        float_value.save()
        models.AnswerKey.objects.create(
            challenge=config.challenge,
            input_element=element,
            value_type=molw_type,
            value_object=float_value,
        )
        return element

    return config, add_element


@pytest.fixture
def benzene_from_mol(molfile_molw_config, elem_factory, float_answer_key_factory, db):
    elem = elem_factory(
        molfile_molw_config.challenge,
        molfile_molw_config.input_type,
        "benzene",
        "ChEBI_16716.mdl",
    )
    answer_key = float_answer_key_factory(
        molfile_molw_config.challenge, elem, molfile_molw_config.output_type, 72.0
    )
    return elem


@pytest.fixture
def smiles_type(challenge, db):
    return models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=True,
        content_type=ContentType.objects.get_for_model(models.TextValue),
        key="smiles",
        description="SMILES",
    )


@pytest.fixture
def molw_type(challenge, db):
    return models.ValueType.objects.create(
        challenge=challenge,
        is_input_flag=False,
        content_type=ContentType.objects.get_for_model(models.FloatValue),
        key="molWeight",
        description="Molecular Weight",
    )


@pytest.fixture
def input_elements(smiles_molw_config, db):
    challenge = smiles_molw_config.challenge
    smiles_type = smiles_molw_config.input_type
    molw_type = smiles_molw_config.output_type
    elems = []
    for idx, (name, smiles) in enumerate(
        [
            ("benzene", "c1ccccc1"),
            ("octane", "C" * 8),
            ("hexane", "C" * 6),
            ("phenol", "c1ccccc1O"),
        ]
    ):
        elem = models.InputElement.objects.create(
            name=name, challenge=challenge, is_public=idx % 2
        )
        smiles_value = models.TextValue.from_string(
            smiles, challenge=challenge, input_element=elem
        )
        smiles_value.save()
        models.InputValue.objects.create(
            input_element=elem, value_type=smiles_type, value_object=smiles_value
        )
        float_value = models.FloatValue.from_string(
            "72.0", challenge=challenge, input_element=elem
        )
        float_value.save()
        models.AnswerKey.objects.create(
            challenge=challenge,
            input_element=elem,
            value_type=molw_type,
            value_object=float_value,
        )
        elems.append(elem)
    return elems


@pytest.fixture
def custom_string_arg(container_arg_factory, draft_submission):
    return container_arg_factory(
        draft_submission.container, key="stringarg", string_value="hello world"
    )


@pytest.fixture
def custom_file_arg(container_arg_factory, draft_submission):
    return container_arg_factory(
        draft_submission.container,
        key="filearg",
        file_name="example.txt",
        file_body="these are the contents of the txt file",
    )
