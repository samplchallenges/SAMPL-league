import csv
import os.path
from collections import defaultdict
from functools import partial

from django.core.files import File
from django.db import models, transaction

from .. import batching, filecache, values_helper
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
        create_file_func = partial(BatchFile.from_local, batch=self)
        create_akeyfile_func = partial(AnswerKeyBatchFile.from_local, input_batch=self)
        with transaction.atomic():
            self._set_elements(input_elements)
            input_value_types = self.batch_group.challenge.valuetype_set.filter(
                is_input_flag=True, on_parent_flag=False
            )
            batching.batchup_elements(
                input_value_types,
                list(self.elements()),
                create_file_func,
                values_helper.get_values,
            )
            output_value_types = self.batch_group.challenge.valuetype_set.filter(
                is_input_flag=False, on_parent_flag=False
            )
            batching.batchup_elements(
                output_value_types,
                list(self.elements()),
                create_akeyfile_func,
                values_helper.answerkey_values,
            )

    def save_score(self, submission_run, score_type, value):
        scores = {}
        with open(value, encoding="utf8") as fp:
            reader = csv.DictReader(fp)
            for row in reader:
                scores[int(row["id"])] = row["value"]
        try:
            for element in self.elements():
                element.save_score(submission_run, score_type, scores[element.id])
        except KeyError as kexc:
            raise ValueError("Missing a score") from kexc

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
    challenge_path = os.path.join("batches", f"{challenge.id}")
    group_path = os.path.join("batch_group", f"{batch_group.id}")
    batch_path = os.path.join("b", f"{instance.batch.id}")
    return os.path.join(
        challenge_path, group_path, batch_path, instance.value_type.key, filename
    )


class BatchFileBase(Timestamped):
    value_type = models.ForeignKey(ValueType, on_delete=models.CASCADE)

    class Meta:
        abstract = True

    @classmethod
    def files_by_key(cls, **filters):
        by_key = {}
        for instance in cls.objects.filter(**filters):
            by_key[instance.value_type.key] = filecache.ensure_local_copy(instance.data)
        return by_key

    @classmethod
    def from_local(cls, filepath, *, value_type, **kwargs):
        filename = os.path.basename(filepath)
        instance = cls(data=filename, value_type=value_type, **kwargs)
        with open(filepath, "rb") as fp:
            instance.data.save(filename, File(fp))
        filecache.preserve_local_copy(instance.data, filepath)
        instance.save()
        return instance


class BatchFile(BatchFileBase):
    batch = models.ForeignKey(InputBatch, on_delete=models.CASCADE)
    data = models.FileField(upload_to=batch_upload_location)

    def __str__(self):
        return f"{self.batch} {self.value_type} {self.data}"


def answer_key_batch_upload_location(instance, filename):
    input_batch = instance.input_batch

    challenge = input_batch.batch_group.challenge
    challenge_path = os.path.join("batch_answer_keys", f"{challenge.id}")
    group_path = os.path.join("batch_group", f"{input_batch.batch_group.id}")
    batch_path = os.path.join("b", f"{input_batch.id}")
    return os.path.join(
        challenge_path, group_path, batch_path, instance.value_type.key, filename
    )


class AnswerKeyBatchFile(BatchFileBase):
    input_batch = models.ForeignKey(InputBatch, on_delete=models.CASCADE)
    data = models.FileField(upload_to=answer_key_batch_upload_location)

    def __str__(self):
        return f"{self.input_batch} {self.value_type} {self.data}"
