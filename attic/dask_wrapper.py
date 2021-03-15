from dask.distributed import Client


class MyDaskClient():
    def __init__(self, address=None):
        self._client = Client(address)

    def _who_has(self, key):
        who_has_dict = self._client.who_has()
        if key in who_has_dict:
            return {"key": key, "worker": who_has_dict[key]}

    def get_status(self, key):
        # first we check if a worker has it
        processing_dict = self._client.processing()
        for worker in processing_dict.keys():
            if key in processing_dict[worker]:
                return {"status": "running", "worker": worker}
        # then we check if the task is in the stream
        for task in reversed(self._client.get_task_stream()):
            if task["key"] == key:
                return {"status": "done", "dask_status": task["status"]}
