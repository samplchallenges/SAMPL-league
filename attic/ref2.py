import logging
import os
from pathlib import Path

from dask.distributed import Client

LOGLEVEL = os.environ.get("LOGLEVEL", "WARNING").upper()
logging.basicConfig(level=LOGLEVEL)


class Submission:
    def __init__(self, challenge_id: int):
        self.challenge_id = challenge_id
        self.challenge_data = CHALLENGE_DB[self.challenge_id]
        self.args = self.challenge_data["args"]
        self.image = self.challenge_data["image"]
        self.scoring_image = self.challenge_data["scoreing_image"]
        self.private_data = self.challenge_data["PRIVATE_DATA"]
        self.public_data = self.challenge_data["PUBLIC_DATA"]
        self.submission_dir = self.challenge_data["SUBMISSION_DIR"]
        self.dask_url = self.challenge_data["DASK_URL"]

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

    @staticmethod
    def _get_mol_list(mol_path):
        if mol_path.exists() and mol_path.is_dir():
            mol_path = mol_path.resolve()
            mols = list(mol_path.glob("mol_*"))
            # return or set self.mols = mols ?
            return mols
        else:
            raise OSError("does not exist or not a directory")

    def _prep_outputs(self):
        run_id = str(self.image_hash)
        # slice off sha:
        run_id = run_id.split(":")[-1]

        # would be nice to have
        # name =  dir.joinpath("foo").mkdir() return path

        self.output_dir_public = self.submission_dir.joinpath(run_id, "public")
        self.output_dir_public.mkdir(parents=True, exist_ok=True)

        self.output_dir_private = self.submission_dir.joinpath(run_id, "private")
        self.output_dir_private.mkdir(parents=True, exist_ok=True)

    def test_container(self):
        mol_list = self._get_mol_list(self.public_data)
        dask_client = self._setup_dask()
        submissions = []
        for mol in mol_list:
            future = dask_client.submit(
                run_dask,
                self.image,
                mol,
                self.output_dir_public,
                self.args,
                use_cache=True,
                pure=False,
            )
            submissions.append(future)
        # We really only need to return/keep track of the keys
        # the code below sums logP calcs -- which means nothing
        print(submissions)
        results = dask_client.gather(submissions)
        print(results)
        total = sum(results)
        print(f"total RMSE for {self.image} is {total}")

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

    def _setup_dask(self):
        return Client(self.dask_url)


def run_dask(image, input_mol, output_dir, args, use_cache=True):
    if use_cache is False:
        raise NotImplementedError
    mol_name = input_mol.name
    try:
        with open(output_dir.joinpath(f"score_{mol_name}"), "r") as score_file:
            score = score_file.readline().strip()
        logging.debug(f"cache hit {mol_name}")
        return float(score)
    except (FileNotFoundError):
        with open(input_mol, "r") as mol_file:
            mol = mol_file.readline().strip()

        command = [mol] + args
        logging.debug(f"running command: {command} on {mol_name} {mol}")
        import docker

        client = docker.from_env()
        result = client.containers.run(image, command, auto_remove=True)
        # Will want to catch an excpetion here and pass the log if done in public
        result = float(result.strip())
        with open(output_dir.joinpath(f"score_{mol_name}"), "w") as score_file:
            score_file.write(f"{str(result)}\n")
        return result


PUBLIC_DATA = Path("data-public").resolve()
PRIVATE_DATA = Path("data-private").resolve()
SUBMISSION_DIR = Path("data-submissions").resolve()

# Mock things that we should pull from DB as dictionary
CHALLENGE_DB = {
    0: {
        "PUBLIC_DATA": PUBLIC_DATA,
        "PRIVATE_DATA": PRIVATE_DATA,
        "SUBMISSION_DIR": SUBMISSION_DIR,
        "DASK_URL": "127.0.0.1:8786",
        "image": "localhost:5000/mmh42/sampl-test:0.1",
        "scoreing_image": "mmh42/score-submission:0.1",
        "challenge_name": "LogP",
        "args": ["--fuzz"],
    }
}

challenge_id = 0
submission = Submission(challenge_id)
submission.prep()
submission.test_container()
# submission.run()
# submission.grade()
