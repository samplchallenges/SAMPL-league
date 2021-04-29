from django.views.generic.detail import DetailView
from django.views.generic.list import ListView

from ..models import Challenge


class ChallengeDetail(DetailView):
    model = Challenge


class ChallengeList(ListView):
    model = Challenge
