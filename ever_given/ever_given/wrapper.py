import docker


def run_container(container_uri, command):
    client = docker.from_env()
    result = client.containers.run(
        container_uri,
        command,
        network_disabled=True,
        network_mode="none",
        auto_remove=True,
    )
    return result
