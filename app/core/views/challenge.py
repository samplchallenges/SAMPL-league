import os.path

import ever_given.wrapper
from django.utils.safestring import mark_safe
from django.views.generic.detail import DetailView
from django.views.generic.list import ListView
from django.utils import timezone

from ..models import Challenge, FileValue


def _input_elements(challenge):
    elements = (
        challenge.inputelement_set.filter(is_public=True)
        .prefetch_related("inputvalue_set", "inputvalue_set__value_type")
        .all()
    )
    input_types = (
        challenge.valuetype_set.filter(is_input_flag=True).order_by("key").all()
    )
    element_attributes = []

    for element in elements:
        input_values = {
            input_value.value_type.key: input_value
            for input_value in element.inputvalue_set.all()
        }
        element_attribute = [element.name]

        for input_type in input_types:
            input_value = input_values[input_type.key]
            if input_type.content_type.model_class() == FileValue:
                value = mark_safe(
                    f'<a href="/download_input/{input_value.pk}/">{input_value.value}</a>'
                )
            else:
                value = input_value.value
        element_attribute.append(value)

        kwargs, file_kwargs = element.all_values()
        args_dict = kwargs
        for k, v in file_kwargs.items():
            args_dict[k] = os.path.basename(v)
        element_attribute.append(ever_given.wrapper.prepare_commandline("", args_dict))
        element_attributes.append(element_attribute)
    return element_attributes, input_types


class ChallengeDetail(DetailView):
    model = Challenge

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        challenge = context["challenge"]
        context["submissions"] = challenge.submission_set.filter(
            user=self.request.user
        ).all()
        context["output_types"] = challenge.valuetype_set.filter(
            is_input_flag=False
        ).order_by("key")
        context["elements"], context["input_types"] = _input_elements(challenge)
        return context


class ChallengeList(ListView):
    model = Challenge
