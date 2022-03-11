class CancelledException(Exception):
    """
    Raised when cancel_requested_func returns True
    """


class LogHandlerBase:
    """
    Subclass this to do something more sophisticated with stdout / stderr
    """

    def handle_stdout(self, log: str) -> None:
        print(log, end="")

    def handle_stderr(self, log: str) -> None:
        print(log, end="")

