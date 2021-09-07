from django.views.generic.detail import DetailView
from django.views.generic.list import ListView

from ..models import Challenge


class ChallengeDetail(DetailView):
    model = Challenge

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["submissions"] = context["challenge"].submission_set.filter(user=self.request.user).all()
        return context

class ChallengeList(ListView):
    model = Challenge
