from django.conf import settings
from django.contrib.auth import authenticate, login, logout
from django.http import HttpResponseBadRequest, HttpResponseForbidden
from django.shortcuts import redirect, render, reverse
from django.views.generic import TemplateView

from ..forms import RegisterForm


class ProfileView(TemplateView):
    template_name = "profile.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)


def register(request):
    if request.method == "GET":
        return render(request, "registration/add.html", {"form": RegisterForm})

    if not settings.ENABLE_REGISTRATION:
        raise HttpResponseForbidden
    if request.method == "POST":
        form = RegisterForm(request.POST)
        if form.is_valid():
            user = form.save()
            login(request, user)
            return redirect(reverse("root"))
        return render(request, "registration/add.html", {"form": RegisterForm})
