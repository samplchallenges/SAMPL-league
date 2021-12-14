"""
Run subprocess with stdout / stderr handling, canceling
"""
import io
import queue
import subprocess
import threading


QUEUE_WAIT_SECONDS = 2
AFTER_TERMINATE_WAIT = 5


class CancelledException(Exception):
    """
    Raised when cancel_requested_func returns True
    """


class LogHandler:
    """
    Subclass this to do something more sophisticated with stdout / stderr
    """

    def handle_stdout(self, log):
        print("stdout", log)

    def handle_stderr(self, log):
        print("stderr", log)


def run(args: list[str], *, log_handler: LogHandler = None, cancel_requested_func=None):
    """
    Runs the command line in args via subprocess.
    Passes stderr and stdout to log_handler and checks whether
    """
    output_buffer = io.BytesIO()
    err_message_queue: queue.Queue = queue.Queue()
    out_message_queue: queue.Queue = queue.Queue()
    if log_handler is None:
        log_handler = LogHandler()

    process = subprocess.Popen(args, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

    out_thread = threading.Thread(
        target=_read_pipe,
        name="stdout",
        args=(process.stdout, out_message_queue, output_buffer),
    )
    err_thread = threading.Thread(
        target=_read_pipe,
        name="stderr",
        args=(process.stderr, err_message_queue),
    )
    err_thread.start()
    out_thread.start()

    while (err_thread.is_alive() or out_thread.is_alive()) or not (
        err_message_queue.empty() and out_message_queue.empty()
    ):
        _handle(log_handler, err_message_queue, out_message_queue)
        if cancel_requested_func is not None:
            _check_cancel(process, cancel_requested_func, log_handler)

    returncode = process.poll()
    if returncode is None:
        raise subprocess.SubprocessError("Subprocess is still running!")
    return output_buffer.getvalue()


def _read_pipe(pipe, msg_queue, output_buffer=None):
    """
    Runs in a thread. Sit on an output pipe (stderr or stdout)
    and throw the messages into a queue for the main thread.
    If output_buffer is set, append them on there too.
    """
    for log in pipe.readline:
        msg_queue.put(log)
        if output_buffer is not None:
            output_buffer.write(log)


def _handle(
    log_handler: LogHandler,
    err_message_queue: queue.Queue,
    out_message_queue: queue.Queue,
):
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


def _check_cancel(process, cancel_requested_func, log_handler):
    should_cancel = cancel_requested_func()
    if should_cancel:
        log_handler.handle_stderr(b"Cancel requested\n")
        process.terminate()
        try:
            process.wait(AFTER_TERMINATE_WAIT)
        except subprocess.TimeoutExpired:
            process.kill()
        raise CancelledException()
