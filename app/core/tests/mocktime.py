from datetime import datetime, timezone

# start_at time is 2020, 1, 1, hour=1, tzinfo=timezone.utc
# end_at time is 2020, 9, 1, hour=1, tzinfo=timezone.utc


def active():
    return datetime(2020, 7, 1, hour=1, tzinfo=timezone.utc)


def inactive_before():
    return datetime(2019, 1, 1, hour=1, tzinfo=timezone.utc)


def inactive_after():
    return datetime(2020, 9, 2, hour=1, tzinfo=timezone.utc)
