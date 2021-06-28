import os.path

import click
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem

CONFORMATION_KEY = "conformation"
MOLW_KEY = "molWeight"
ATOMCOUNT_KEY = "numAtoms"
BONDCOUNT_KEY = "numBonds"

@click.group()
def cli():
    pass

@cli.command()
@click.option("--conformation_prediction")
@click.option("--conformation_answerkey")
def score_evaluation(conformation_prediction, conformation_answerkey):
    # Actually, load the files, calc RMSE, etc
    print("RMSE 1.345")

    #mol = Chem.MolFromMolFile(molfile)

@cli.command()
@click.option("--RMSE")
def score_submissionrun(rmse):
    pass
