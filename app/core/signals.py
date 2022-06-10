from django.core.exceptions import ValidationError
from django.db.models.signals import post_save, pre_save
from django.dispatch import receiver

from . import batcher, models


@receiver(post_save, sender=models.Challenge)
def manage_batch_group(
    sender, instance, created, raw, using, update_fields, **kwargs
):  # pylint: disable=unused-argument
    if created:
        # No batch group needed when the challenge is new
        return
    if update_fields is not None:
        if "max_batch_size" not in update_fields:
            return
    if instance.max_batch_size == 0:
        return
    # Also verify that max_batch_size is changed
    batcher.generate_batches(instance)


@receiver(pre_save, sender=models.InputBatchMembership)
def _validate_batch_membership(
    sender, instance, **kwargs
):  # pylint: disable=unused-argument
    if instance.batch.batch_group != instance.batch_group:
        raise ValidationError(
            "batch group on InputBatchMembership object doesn't match the"
            "batch group of its batch"
        )
    if instance.batch.is_public != instance.input_element.is_public:
        raise ValidationError(
            f"Input element is_public ({instance.input_element.is_public}) doesn't"
            f" match batch is_public ({instance.batch.is_public})"
        )
