from django.urls import path, include

from .views.challenge import ChallengeDetail, ChallengeList

from .views.submission import (
    SubmissionDetail,
    SubmissionList,
    SubmissionDelete,
    edit_submission_view,
)

from .views.root import IndexView
from .views.profile import ProfileView, register


urlpatterns = [
    path("", IndexView.as_view(), name="root"),
    path("accounts/", include("django.contrib.auth.urls")),
    path("profile/", ProfileView.as_view(), name="profile-view"),
    path("profile/add", register, name="profile-register"),
    path("challenge/", ChallengeList.as_view(), name="challenge-list"),
    path("challenge/<int:pk>/", ChallengeDetail.as_view(), name="challenge-detail"),
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
