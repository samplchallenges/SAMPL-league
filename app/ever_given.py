def run_submission_container(container_uri, command):
    import docker
    client = docker.from_env()
    result = client.containers.run(container_uri, command, auto_remove=True)
    return result
