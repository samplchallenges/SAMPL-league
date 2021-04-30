import django

django.setup()

from referee.tasks import run_submission, score_submission


run_submission(19, is_public=True)
run_submission(19, is_public=False)

score_submission(19)
# fire_off_tasks(5, ["c1ccccc1", "C", "CC1=CC=CC=C1"])
