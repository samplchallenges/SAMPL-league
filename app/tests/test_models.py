import pytest
from django.core.exceptions import ValidationError
from django.db import models

from core.models import GenericValue, register_value_model


def test_value_registration():
    with pytest.raises(ValidationError):

        @register_value_model
        class InvalidValueModel(GenericValue):
            class Meta:
                app_label = "core.apps.CoreConfig"

    @register_value_model
    class CharValueModel(GenericValue):
        value = models.CharField(max_length=100)

        class Meta:
            app_label = "core.apps.CoreConfig"
