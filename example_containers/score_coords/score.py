import json
import os.path
import random
import statistics

import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem


@click.group()
def cli():
    pass


@cli.command()
@click.option("--conformation_prediction")
@click.option("--conformation_answerkey")
def score_evaluation(conformation_prediction, conformation_answerkey):
    # Actually, load the files, calc RMSE, etc
    rmse = 10.0 + random.random()*5
    print("RMSE", rmse)

    #mol = Chem.MolFromMolFile(molfile)


@cli.command()
@click.option("--scores")
def score_submissionrun(scores_file):
    with open(scores_file, "r") as fp:
        scores_dicts = json.load(fp)

    rmses = [score_dict["RMSE"] for score_dict in scores_dict]

    print("max_RMSE", max(rmses))
    print("min_RMSE", min(rmses))
    print("median_RMSE", statistics.median(rmses))
    print("mean_RMSE", statistics.mean(rmses))
