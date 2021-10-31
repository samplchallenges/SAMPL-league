import os.path

from django import template
from django.utils.html import format_html

from core.models import FileValue, InputValue

register = template.Library()


@register.filter
def render_prediction(value_parent):
    value = value_parent.value
    if isinstance(value_parent.value_object, FileValue):
        value = format_html(
            '<a href="/download_output/{}/">{}</a>', value_parent.pk, os.path.basename(value.name)
        )
    return value


@register.filter
def render_input_value(input_value):
    if not isinstance(input_value, InputValue):
        return input_value
    value = input_value.value
    if isinstance(input_value.value_object, FileValue):
        value = format_html(
            '<a href="/download_input/{}/">{}</a>',
            input_value.pk, os.path.basename(value.name))
    return value
