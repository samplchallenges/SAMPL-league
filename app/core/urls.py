from django.urls import path

from .views.submission import (
    SubmissionDetail,
    SubmissionList,
    SubmissionDelete,
    edit_submission_view,
)

from .views.root import IndexView
from .views.auth import logout_view


urlpatterns = [
    path("", IndexView.as_view(), name="root"),
    path("logout/", logout_view, name="logout"),
    path("submission/", SubmissionList.as_view(), name="submission-list"),
    path("submission/add", edit_submission_view, name="submission-add"),
    path("submission/<int:pk>/", SubmissionDetail.as_view(), name="submission-detail"),
    path(
        "submission/<int:pk>/edit/",
        edit_submission_view,
        name="submission-update",
    ),
    path(
        "submission/<int:pk>/delete/",
        SubmissionDelete.as_view(),
        name="submission-delete",
    ),
]
