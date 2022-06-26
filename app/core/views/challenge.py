import os.path
from collections import OrderedDict, defaultdict

import ever_given.wrapper
from django.views.generic.detail import DetailView
from django.views.generic.list import ListView

from .. import values_helper
from ..models import Challenge, InputValue


def _input_elements(challenge):
    input_values = (
        InputValue.objects.select_related(
            "input_element", "value_type", "value_type__content_type"
        )
        .filter(
            input_element__challenge=challenge,
            value_type__is_input_flag=True,
            input_element__is_public=True,
        )
        .order_by("input_element_id", "value_type__key")
        .all()
    )
    input_types = OrderedDict()
    elements = defaultdict(OrderedDict)
    for input_value in input_values:
        key = input_value.value_type.key
        if key not in input_types:
            input_types[key] = input_value.value_type
        element_id = input_value.input_element_id
        elements[element_id]["name"] = input_value.input_element.name
        elements[element_id][key] = input_value

    for element_dict in elements.values():
        element_values = [v for k, v in element_dict.items() if k != "name"]
        (
            kwargs,
            file_kwargs,
        ) = values_helper._element_values(  # pylint: disable=protected-access
            element_values
        )
        args_dict = kwargs
        for k, v in file_kwargs.items():
            args_dict[k] = os.path.basename(v)
        element_dict["container_args"] = " ".join(
            ever_given.wrapper.prepare_command_list("", args_dict)
        )
    return list(elements.values()), list(input_types.values())


class ChallengeDetail(DetailView):
    model = Challenge

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        challenge = context["challenge"]
        if self.request.user.is_authenticated:
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
