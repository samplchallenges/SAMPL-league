import docker


def run_container(container_uri, command, inputdir=None, outputdir=None):
    client = docker.from_env()
    volumes = {}
    if inputdir:
        volumes[inputdir] = {
            'bind': '/mnt/inputs',
            'mode': 'ro'}
    if outputdir:
        volumes[outputdir] = {
            'bind': '/mnt/outputs',
            'mode': 'rw'}

    result = client.containers.run(
        container_uri,
        command,
        network_disabled=True,
        network_mode="none",
        #auto_remove=True, TODO: figure out how to get logging back
        # while keeping auto_remove
        volumes=volumes
    )
    return result
