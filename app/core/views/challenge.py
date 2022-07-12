from django.views.generic.detail import DetailView
from django.views.generic.list import ListView

from .. import template_helpers
from ..models import Challenge


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
        context["element_table"] = template_helpers.ElementTable(challenge)
        return context


class ChallengeList(ListView):
    model = Challenge
