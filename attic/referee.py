# First we want to run the submission on the public data + report errors (per mol)
# Run the score container, report results
# Then we run the submission on the private data (per mol)
# score container


def gather_mols():
    pass

def pull_container():
    # We will want to make sure we can the pull the container and can fail fast
    pass

def run_submission(image, mol_path, args):
    import docker
    client = docker.from_env()
    with open(mol_path) as mol_file:
        mol = mol_file.readline().strip()

    command = [mol] + args
    print(f"running command: {command} on mol {mol}")
    result = client.containers.run(image, command, auto_remove=True)
    return result

image = "mmh42/sampl-test:0.1"
mol_path = "data-public/mol_1"
args = ["--fuzz"]
result = run_submission(image, mol_path, args)
print(result)
