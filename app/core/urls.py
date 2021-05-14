from django.urls import include, path

from .views.challenge import ChallengeDetail, ChallengeList
from .views.profile import ProfileView, register
from .views.root import IndexView
from .views.submission import (
    SubmissionDelete,
    SubmissionDetail,
    SubmissionList,
    clone_submission_view,
    edit_submission_view,
    submit_submission_view,
)

urlpatterns = [
    path("", IndexView.as_view(), name="root"),
    path("accounts/", include("django.contrib.auth.urls")),
    path("profile/", ProfileView.as_view(), name="profile-view"),
    path("profile/add/", register, name="profile-register"),
    path("challenge/", ChallengeList.as_view(), name="challenge-list"),
    path("challenge/<int:pk>/", ChallengeDetail.as_view(), name="challenge-detail"),
    path("submission/", SubmissionList.as_view(), name="submission-list"),
    path("submission/<int:pk>/", SubmissionDetail.as_view(), name="submission-detail"),
    path(
        "submission/<int:pk>/submit", submit_submission_view, name="submission-submit"
    ),
    path("submission/add/", edit_submission_view, name="submission-add"),
    path(
        "submission/<int:pk>/clone/",
        edit_submission_view,
        name="submission-clone",
        kwargs={"clone": True},
    ),
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
