import docker


def run_container(container_uri, command, inputdir=None, outputdir=None):
    client = docker.from_env()
    volumes = {}
    if inputdir:
        volumes[inputdir] = {
            'bind': '/mnt/inputs',
            'mode': 'ro'}
    if outputdir:
        outputdir_bind = "/mnt/outputs"
        volumes[outputdir] = {
            'bind': outputdir_bind,
            'mode': 'rw'}
        command = f" --output-dir {outputdir_bind} {command}"

    result = client.containers.run(
        container_uri,
        command,
        volumes=volumes,
        network_disabled=True,
        network_mode="none",
        remove=True,
    )
    return result
