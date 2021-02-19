from pathlib import Path
from typing import List
import dask

# First we want to run the submission on the public data
# + report errors (per mol)
# Run the score container, report results
# Then we run the submission on the private data (per mol)
# score container

SUBMISSION_DIR = Path("data-submissions").resolve()

print(SUBMISSION_DIR)


def gather_mols(data_dir: str) -> List[str]:
    p = Path(data_dir)
    if p.exists() and p.is_dir():
        p = p.resolve()
        mols = list(p.glob("mol_*"))
        return mols
    else:
        raise OSError("does not exist or not a directory")


def prep_out_path(image):
    run_id = get_container_id(image)
    run_id = str(run_id)
    # slice off sha:
    run_id = run_id.split(":")[-1]
    SUBMISSION_DIR.joinpath(run_id).mkdir(parents=True, exist_ok=True)
    return run_id


def get_container_id(image):
    import docker
    client = docker.from_env()
    container = client.images.pull("localhost:5000/mmh42/sampl-test", tag="0.1")
    return container.id

    # We will want to make sure we can the pull the container and can fail fast


@dask.delayed
def run_submission(image: str, mol_path: str, args: List[str], run_id: str) -> float:
    mol_name = mol_path.name
    try:
        with open(SUBMISSION_DIR.joinpath(run_id, f"score_{mol_name}"), "r") as score_file:
            score = score_file.readline().strip()
        print(f"cache hit {mol_name}")
        return float(score)
    except(FileNotFoundError):
        import docker
        client = docker.from_env()
        with open(mol_path, "r") as mol_file:
            mol = mol_file.readline().strip()

        command = [mol] + args
        print(f"running command: {command} on {mol_name} {mol}")
        result = client.containers.run(image, command, auto_remove=True)
        result = float(result.strip())
        with open(SUBMISSION_DIR.joinpath(run_id, f"score_{mol_name}"), "w") as score_file:
            score_file.write(f"{str(result)}\n")
        return result


image = "mmh42/sampl-test:0.1"
args = ["--fuzz"]
mols = gather_mols("data-public")

container_id = get_container_id(image)
print(container_id)
scores = []
run_id = prep_out_path(image)
for mol_path in mols:
    result = run_submission(image, mol_path, args, run_id)
    scores.append(result)
total = dask.delayed(sum)(scores)
total.visualize()
print("computing")
rmse_sum = total.compute()
print(f"total RMSE for {image} is {rmse_sum}")
