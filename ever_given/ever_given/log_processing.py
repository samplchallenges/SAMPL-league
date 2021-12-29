import io
import queue
import threading
import time
import typing

QUEUE_WAIT_SECONDS = 2

from .engines.utils import ContainerInstance
from .utils import CancelledException, LogHandlerBase


def process_messages(
    running_container: ContainerInstance,
    log_handler: LogHandlerBase,
    cancel_requested_func: typing.Optional[typing.Callable[[], bool]],
) -> str:

    output_buffer = io.StringIO()
    err_message_queue: "queue.Queue[str]" = queue.Queue()
    out_message_queue: "queue.Queue[str]" = queue.Queue()
    out_thread = threading.Thread(
        target=_read_stdout,
        name="stdout",
        args=(running_container, out_message_queue, output_buffer),
    )
    err_thread = threading.Thread(
        target=_read_stderr, name="stderr", args=(running_container, err_message_queue)
    )
    err_thread.start()
    out_thread.start()

    while (err_thread.is_alive() or out_thread.is_alive()) or not (
        err_message_queue.empty() and out_message_queue.empty()
    ):
        _handle(log_handler, err_message_queue, out_message_queue)
        if cancel_requested_func is not None:
            _check_cancel(running_container, cancel_requested_func, log_handler)

    return output_buffer.getvalue()


def _read_stdout(
    container: ContainerInstance,
    out_message_queue: queue.Queue,
    output_buffer: io.StringIO,
):
    for log in container.logs(want_stdout=True, want_stderr=False):
        out_message_queue.put(log)
        output_buffer.write(log)


def _read_stderr(container: ContainerInstance, err_message_queue: queue.Queue):
    for log in container.logs(want_stdout=False, want_stderr=True):
        err_message_queue.put(log)


def _handle(
    log_handler: LogHandlerBase,
    err_message_queue: queue.Queue,
    out_message_queue: queue.Queue,
) -> None:
    """
    Runs from a loop in the main thread. Look in message queues for inputs
    and pass them to log_handler
    Timeout after a couple of seconds on each queue so we  can check the other queue
    and then check whether we should cancel.
    """
    try:
        while True:
            log = out_message_queue.get(timeout=QUEUE_WAIT_SECONDS)
            if log:
                log_handler.handle_stdout(log)
    except queue.Empty:
        pass
    try:
        while True:
            log = err_message_queue.get(timeout=QUEUE_WAIT_SECONDS)
            if log:
                log_handler.handle_stderr(log)
    except queue.Empty:
        pass


def _check_cancel(
    container: ContainerInstance,
    cancel_requested_func: typing.Callable[[], bool],
    log_handler: LogHandlerBase,
) -> None:
    should_cancel = cancel_requested_func()
    if should_cancel:
        log_handler.handle_stderr("Cancel requested\n")
        container.kill()
        raise CancelledException()
