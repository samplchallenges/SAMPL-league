import time

from ever_given import wrapper

container_uri = "robbason/logging-example:latest"

final_command = wrapper._prepare_commandline("", {"smiles": "c1ccccc1"})
results = wrapper.run_container(container_uri, final_command, {}, output_dir="")

# container_uri, file_kwargs={}, kwargs={"smiles": "c1ccccc1"})
for x in range(45):
    print("getting logs")
    print("stderr", results.logs(stdout=False))
    print("stdout", results.logs(stderr=False))
    results.reload()
    print("status", results.status)
    if results.status == "exited":
        break
    time.sleep(1)

print("ALL DONE")
