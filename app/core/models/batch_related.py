import os.path
import tempfile
from collections import defaultdict

from django.core.files import File
from django.db import models, transaction

from .. import filecache
from ..batching import BATCHERS
from .admin_managed import Challenge, InputElement, ValueType
from .infra_models import Timestamped


class InputBatchGroup(Timestamped):
    """
    When batch sizes are changed, we want to keep history associated with old batches
    """

    challenge = models.ForeignKey(Challenge, on_delete=models.CASCADE)
    max_batch_size = models.PositiveIntegerField(
        help_text="Value of setting when batch group was created"
    )

    def generate_batches(self):
        for is_public in (True, False):
            elements = self.challenge.inputelement_set.filter(
                is_public=is_public, is_parent=False
            ).order_by("parent_id", "name")
            unparented = []
            by_parent = defaultdict(list)
            for element in elements:
                if element.parent_id is None:
                    unparented.append(element)
                else:
                    by_parent[element.parent_id].append(element)
            for parent_id, elements in ((None, unparented), *by_parent.items()):
                if elements:
                    for idx in range(0, len(elements), self.max_batch_size):
                        batch_elements = elements[idx : idx + self.max_batch_size]
                        batch = self.inputbatch_set.create(
                            is_public=is_public, parent_input_element_id=parent_id
                        )
                        batch.batchup(batch_elements)

    def __str__(self):
        return f'Challenge "{self.challenge.name}" batch size {self.max_batch_size}'


class InputBatch(Timestamped):
    batch_group = models.ForeignKey(InputBatchGroup, on_delete=models.CASCADE)
    parent_input_element = models.ForeignKey(
        InputElement, on_delete=models.CASCADE, null=True, blank=True
    )
    is_public = models.BooleanField(default=False)

    def _set_elements(self, input_elements):
        self.inputbatchmembership_set.all().delete()
        for input_element in input_elements:
            if self.parent_input_element_id != input_element.parent_id:
                raise ValueError(
                    "All input elements to a batch must share a "
                    "parent id (or all have no parent elements)"
                )
            InputBatchMembership.objects.create(
                batch=self, batch_group=self.batch_group, input_element=input_element
            )

    def elements(self):
        for ibm in self.inputbatchmembership_set.select_related("input_element"):
            yield ibm.input_element

    def batchup(self, input_elements):
        with transaction.atomic():
            self._set_elements(input_elements)
            with tempfile.TemporaryDirectory() as temp_dir:
                for value_type in self.batch_group.challenge.valuetype_set.filter(
                    is_input_flag=True, on_parent_flag=False
                ):
                    batcher = BATCHERS.get(value_type.batch_method)
                    if batcher is None:
                        raise ValueError(
                            f"Batch method must be set on value type {value_type.key}"
                        )
                    output_path = os.path.join(
                        temp_dir, f"{value_type.key}.{batcher.suffix}"
                    )
                    batcher.call(input_elements, value_type.key, output_path)
                    batch_file = BatchFile.from_local(
                        output_path, batch=self, value_type=value_type
                    )
                    batch_file.save()

    def __str__(self):
        return f"Public? {self.is_public}"


class InputBatchMembership(Timestamped):
    batch_group = models.ForeignKey(InputBatchGroup, on_delete=models.CASCADE)
    batch = models.ForeignKey(InputBatch, on_delete=models.CASCADE)
    input_element = models.ForeignKey(InputElement, on_delete=models.CASCADE)

    class Meta:
        unique_together = ["batch_group", "input_element"]

    def __str__(self):
        return str(self.input_element)


def batch_upload_location(instance, filename):
    batch_group = instance.batch.batch_group
    challenge = batch_group.challenge
    group_path = os.path.join("batch_group", f"{batch_group.id}")
    return os.path.join(
        "batches", f"{challenge.id}", group_path, instance.value_type.key, filename
    )


class BatchFile(Timestamped):
    batch = models.ForeignKey(InputBatch, on_delete=models.CASCADE)
    value_type = models.ForeignKey(ValueType, on_delete=models.CASCADE)
    data = models.FileField(upload_to=batch_upload_location)

    def __str__(self):
        return f"{self.batch} {self.value_type} {self.data}"

    @classmethod
    def from_local(cls, filepath, *, batch, value_type):
        filename = os.path.basename(filepath)
        instance = cls(data=filename, batch=batch, value_type=value_type)
        with open(filepath, "rb") as fp:
            instance.data.save(filename, File(fp))
        filecache.preserve_local_copy(instance.data, filepath)
        return instance
