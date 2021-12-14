import docker

def run_container(container_uri, command, inputdir_map=None, output_dir=None):
    client = docker.from_env()
    volumes = {}
    for inputdir, guest_input_dir in inputdir_map.items():
        volumes[str(inputdir)] = {"bind": str(guest_input_dir), "mode": "ro"}
    if output_dir:
        output_dir = Path(output_dir).resolve()
        volumes[str(output_dir)] = {"bind": str(GUEST_OUTPUT_DIR), "mode": "rw"}
        command = f" --output-dir {GUEST_OUTPUT_DIR} {command}"

    running_container = client.containers.run(
        container_uri,
        command,
        volumes=volumes,
        network_disabled=True,
        network_mode="none",
        remove=False,
        detach=True,
    )
    return running_container
