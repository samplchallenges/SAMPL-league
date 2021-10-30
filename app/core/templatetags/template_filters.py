from django import template
from django.utils.html import format_html

from core.models import FileValue

register = template.Library()


@register.filter
def render_prediction(value_parent):
    value = value_parent.value
    if isinstance(value_parent.value_object, FileValue):
        value = format_html('<a href="/download_output/{}/">{}</a>',
                            value_parent.pk, value)
    return value
