import time

import pytest

import ever_given.wrapper
from ever_given.log_processing import QUEUE_WAIT_SECONDS, CancelledException


def test_cancellation():
    # Use a slow container so we have time to cancel it
    container_uri = "ghcr.io/megosato/logging-example:latest"
    kwargs = {"smiles": "c1cccnc1"}
    start_at = time.time()
    with pytest.raises(CancelledException):
        results = {
            key: value
            for key, value in ever_given.wrapper.run(
                container_uri,
                kwargs=kwargs,
                file_kwargs={},
                cancel_requested_func=lambda: True,
            )
        }
    end_at = time.time()
    # Should cancel in not much more time than QUEUE_WAIT_SECONDS * 2
    # in log_processing
    assert end_at - start_at < QUEUE_WAIT_SECONDS * 2 + 2


def test_no_cancellation():
    container_uri = "ghcr.io/megosato/calc-molwt:latest"
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri,
            kwargs=kwargs,
            file_kwargs={},
            cancel_requested_func=lambda: False,
        )
    }
