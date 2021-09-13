import io
import queue
import threading


def process_messages(running_container, log_handler):
    if log_handler is None:
        log_handler = PrintLogHandler()

    if not (
        hasattr(log_handler, "handle_stdout") and hasattr(log_handler, "handle_stderr")
    ):
        raise TypeError("log_handler must have handle_stdout and handle_stderr methods")

    output_buffer = io.BytesIO()
    err_message_queue = queue.Queue()
    out_message_queue = queue.Queue()
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

    return output_buffer.getvalue()


class PrintLogHandler:
    def handle_stdout(self, log):
        print("stdout", log)

    def handle_stderr(self, log):
        print("stderr", log)


def _read_stdout(container, out_message_queue, output_buffer):
    for log in container.logs(stdout=True, stderr=False, stream=True):
        out_message_queue.put(log)
        output_buffer.write(log)


def _read_stderr(container, err_message_queue):
    for log in container.logs(stdout=False, stderr=True, stream=True):
        err_message_queue.put(log)


def _handle(log_handler, err_message_queue, out_message_queue):
    QUEUE_WAIT_SECONDS = 2

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
