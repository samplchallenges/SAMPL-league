from pathlib import Path
from typing import List
import dask

# First we want to run the submission on the public data
# + report errors (per mol)
# Run the score container, report results
# Then we run the submission on the private data (per mol)
# score container


def gather_mols(data_dir: str) -> List[str]:
    p = Path(data_dir)
    if p.exists() and p.is_dir():
        p = p.resolve()
        mols = list(p.glob("mol_*"))
        return mols
    else:
        raise OSError("does not exist or not a directory")


def get_container_id(image):
    import docker
    client = docker.from_env()
    container = client.images.pull("localhost:5000/mmh42/sampl-test", tag="0.1")
    return container.id

    # We will want to make sure we can the pull the container and can fail fast


@dask.delayed
def run_submission(image: str, mol_path: str, args: List[str]) -> float:
    import docker
    client = docker.from_env()
    mol_name = mol_path.name
    with open(mol_path) as mol_file:
        mol = mol_file.readline().strip()

    command = [mol] + args
    print(f"running command: {command} on {mol_name} {mol}")
    result = client.containers.run(image, command, auto_remove=True)
    result = float(result.strip())
    return result


image = "mmh42/sampl-test:0.1"
args = ["--fuzz"]
mols = gather_mols("data-public")

container_id = get_container_id(image)
print(container_id)
scores = []
for mol_path in mols:
    result = run_submission(image, mol_path, args)
    scores.append(result)
total = dask.delayed(sum)(scores)
total.visualize()
print("computing")
total.compute()
