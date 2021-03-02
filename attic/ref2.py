from pathlib import Path

from dask.distributed import Client


class Submission:
    def __init__(self, challenge_id: int):
        self.challenge_id = challenge_id
        self.challenge_data = CHALLENGE_DB[self.challenge_id]
        self.args = self.challenge_data["args"]
        self.image = self.challenge_data["image"]
        self.scoring_image = self.challenge_data["scoreing_image"]

    def prep(self):
        self._pull_image()
        self._prep_outputs()

    def _pull_image(self):
        # Prime docker to have the image
        import docker

        self.client = docker.from_env()
        # Might need to add a wait OR report progress on this
        # Could be slow for large images
        # Could also error out from http errors, need to retry
        container = self.client.images.pull(self.image)
        self.image_hash = container.attrs.get("RepoDigests")[0]

    def _get_mol_list(self, mol_path):
        if mol_path.exists() and mol_path.is_dir():
            mol_path = mol_path.resolve()
            mols = list(mol_path.glob("mol_*"))
            self.mols = mols
        else:
            raise OSError("does not exist or not a directory")

    def _prep_outputs(self):
        run_id = str(self.image_hash)
        # slice off sha:
        run_id = run_id.split(":")[-1]
        self.output_dir_public = SUBMISSION_DIR.joinpath(run_id, "public").mkdir(
            parents=True, exist_ok=True
        )
        self.output_dir_public = SUBMISSION_DIR.joinpath(run_id, "private").mkdir(
            parents=True, exist_ok=True
        )

    def _run(self, input_mol, output_dir, cache):
        mol_name = input_mol.name
        try:
            with open(output_dir.joinpath(f"score_{mol_name}"), "r") as score_file:
                score = score_file.readline().strip()
            print(f"cache hit {mol_name}")
            return float(score)
        except (FileNotFoundError):
            with open(input_mol, "r") as mol_file:
                mol = mol_file.readline().strip()

            command = [mol] + self.args
            print(f"running command: {command} on {mol_name} {mol}")
            result = self.client.containers.run(self.image, command, auto_remove=True)
            # Will want to catch an excpetion here and pass the log if done in public
            result = float(result.strip())
            with open(output_dir.joinpath(f"score_{mol_name}"), "w") as score_file:
                score_file.write(f"{str(result)}\n")
            return result

    def test_container(self):
        # run on public data + score
        # also return the logs publicly
        pass

    def run(self):
        # run on private data + score
        # use context manager
        pass

    def _score(self):
        pass

    def grade(self):
        # cleanup stuff
        # update db with values?
        pass


Client

PUBLIC_DATA = Path("data-public").resolve()
PRIVATE_DATA = Path("data-private").resolve()
SUBMISSION_DIR = Path("data-submissions").resolve()

# Mock things that we should pull from DB as dictionary
CHALLENGE_DB = {
    0: {
        "image": "localhost:5000/mmh42/sampl-test:0.1",
        "scoreing_image": "mmh42/score-submission:0.1",
        "challenge_name": "LogP",
        "args": ["--fuzz"],
    }
}

challenge_id = 0
submission = Submission(challenge_id)
submission.prep()
# submission.test_container()
# submission.run()
# submission.grade()
