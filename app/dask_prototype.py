import django

django.setup()

from referee.tasks import get_status, run_submission, score_submission

srun = run_submission(18, is_public=True)
status = get_status(srun.digest)
print(status)
# run_submission(19, is_public=False)

# score_submission(19)
# fire_off_tasks(5, ["c1ccccc1", "C", "CC1=CC=CC=C1"])
