from referee.tasks import run_submission, score_submission

run_submission(5, is_public=True)
run_submission(5, is_public=False)

score_submission(5)
# fire_off_tasks(5, ["c1ccccc1", "C", "CC1=CC=CC=C1"])
