from django.urls import path

from .views.submission import (
    SubmissionDetail,
    SubmissionList,
    SubmissionCreate,
    SubmissionUpdate,
    SubmissionDelete,
)


urlpatterns = [
    path("submission/", SubmissionList.as_view(), name="submission-list"),
    path("submission/add", SubmissionCreate.as_view(), name="submission-add"),
    path("submission/<int:pk>/", SubmissionDetail.as_view(), name="submission-detail"),
    path(
        "submission/<int:pk>/edit/",
        SubmissionUpdate.as_view(),
        name="submission-update",
    ),
    path(
        "submission/<int:pk>/delete/",
        SubmissionDelete.as_view(),
        name="submission-delete",
    ),
]
