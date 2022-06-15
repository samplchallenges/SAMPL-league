"""
Given an InputElement, load pairs of dicts of value_type.key: value
one for files, one for non-files
"""
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    # Avoid circular import
    # TODO - in Python 3.11, the type annotations won't need the double-quotes
    # see https://peps.python.org/pep-0563/#abstract
    from .models.values import FileValue, InputElement, InputValue

from . import filecache


def _element_values(input_values: List["InputValue"]):
    file_values = {}
    values = {}

    for input_value in input_values:
        value_type = input_value.value_type
        key = value_type.key
        content_type = value_type.content_type
        if content_type.model_class() == FileValue:
            file_values[key] = filecache.ensure_local_copy(input_value.value)
        else:
            values[key] = input_value.value
    return values, file_values


def _values(input_element: "InputElement"):
    """
    Returns a pair of key: value dicts, where the first dict is regular values
    and the second is file values
    """
    input_values = input_element.inputvalue_set.select_related(
        "value_type", "value_type__content_type"
    ).all()
    return _element_values(input_values)


def all_values(input_element: "InputElement"):
    """
    Includes parent values
    """
    values, file_values = _values(input_element)
    if input_element.parent:
        (
            parent_values,
            parent_file_values,
        ) = input_element.parent._values()  # pylint: disable=protected-access
        values.update(parent_values)
        file_values.update(parent_file_values)
    return values, file_values
