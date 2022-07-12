# pylint: disable=unused-argument, unused-variable, protected-access

import pytest
from django.core.exceptions import ValidationError
from django.db import IntegrityError

from core import models


def test_batch_build_works(smiles_molw_config, input_elements):
    group1 = models.InputBatchGroup.objects.create(
        challenge=smiles_molw_config.challenge, max_batch_size=1
    )
    batch1_priv = models.InputBatch.objects.create(batch_group=group1, is_public=False)
    batch1_pub = models.InputBatch.objects.create(batch_group=group1, is_public=True)

    public_elements = [element for element in input_elements if element.is_public]
    private_elements = [element for element in input_elements if not element.is_public]

    batch1_pub.batchup(public_elements)
    batch1_priv.batchup(private_elements)


def test_generate_batches(smiles_molw_config, input_elements):
    challenge = smiles_molw_config.challenge
    challenge.max_batch_size = 2
    challenge.save()
    batch_group = challenge.current_batch_group()
    assert batch_group.inputbatch_set.count() == 2
    batch = batch_group.inputbatch_set.filter(is_public=True).get()
    assert batch.inputbatchmembership_set.count() == 2
    assert batch.batchfile_set.count() == 1


def test_batch_checks(smiles_molw_config, input_elements):
    element = input_elements[1]
    assert element.is_public
    group = models.InputBatchGroup.objects.create(
        challenge=smiles_molw_config.challenge, max_batch_size=1
    )
    batch_priv = models.InputBatch.objects.create(batch_group=group)
    batch1_pub = models.InputBatch.objects.create(batch_group=group, is_public=True)
    batch2_pub = models.InputBatch.objects.create(batch_group=group, is_public=True)
    with pytest.raises(ValidationError):
        batch_priv._set_elements([element])

    batch1_pub._set_elements([element])

    with pytest.raises(IntegrityError):
        batch2_pub._set_elements([element])
