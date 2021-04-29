import pytest
from django.core.exceptions import ValidationError
from django.db import models

from core.models import GenericOutputValue, Prediction


def test_output_value_registration():
    with pytest.raises(ValidationError):

        @Prediction.register_value_model
        class InvalidValueModel(GenericOutputValue):
            class Meta:
                app_label = "core.apps.CoreConfig"

    @Prediction.register_value_model
    class CharValueModel(GenericOutputValue):
        value = models.CharField(max_length=100)

        class Meta:
            app_label = "core.apps.CoreConfig"
