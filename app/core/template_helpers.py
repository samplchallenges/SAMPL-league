import os.path
from collections import OrderedDict, defaultdict

import ever_given.wrapper

from . import values_helper
from .models import InputValue


class ElementTable:
    child_types = OrderedDict()
    parent_types = OrderedDict()
    elems = defaultdict(OrderedDict)
    parent_elems = defaultdict(OrderedDict)
    rows = []

    def __init__(self, challenge):
        input_values = (
            InputValue.objects.select_related(
                "input_element",
                "value_type",
                "value_type__content_type",
            )
            .filter(
                input_element__challenge=challenge,
                value_type__is_input_flag=True,
                input_element__is_public=True,
            )
            .order_by("input_element_id", "value_type__key")
            .all()
        )

        for input_value in input_values:
            key = input_value.value_type.key
            elem = input_value.input_element
            if elem.is_parent:
                holder = self.parent_elems
                type_holder = self.parent_types
            else:
                holder = self.elems
                type_holder = self.child_types
                holder[elem.id]["parent"] = elem.parent_id

            if key not in type_holder:
                type_holder[key] = input_value.value_type

            holder[elem.id]["name"] = elem.name
            holder[elem.id][key] = input_value
        for elem in self.elems.values():
            self.rows.append(self.__collate_row(elem))

    def __prep_cmdline(self, element_values):
        (
            kwargs,
            file_kwargs,
        ) = values_helper.element_values(element_values)
        args_dict = kwargs
        for k, v in file_kwargs.items():
            args_dict[k] = os.path.basename(v)
        return " ".join(ever_given.wrapper.prepare_command_list("", args_dict))

    def __collate_row(self, elem):
        element_values = [v for k, v in elem.items() if k not in ("name", "parent")]
        parent_id = elem["parent"]
        if parent_id:
            parent_values = [
                v for k, v in self.parent_elems[parent_id].items() if k not in ("name")
            ]
            element_values = parent_values + element_values
        container_args = self.__prep_cmdline(element_values)
        return [elem["name"], *element_values, container_args]

    @property
    def types(self):
        return list(self.parent_types.values()) + list(self.child_types.values())
