import docker


def run_submission_container(container_uri, command):
    client = docker.from_env()
    result = client.containers.run(container_uri, command, auto_remove=True)
    return result
