import json
import os.path
import random
import statistics
import argparse

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem
from rdkit.Chem import AllChem


def score_evaluation(conformation_prediction, conformation_answerkey,
                     molweight_prediction, molweight_answerkey):
    # Actually, load the files, calc RMSE, etc
    rmse = 10.0 + random.random()*5
    print("rmse", rmse)
    print("diff", rmse)

    #mol = Chem.MolFromMolFile(molfile)


def score_submissionrun(scores):
    scores_file = scores
    with open(scores_file, "r") as fp:
        scores_dict = json.load(fp)

    rmses = [score_dict["diff"] for score_dict in scores_dict]

    print("max_RMSE", max(rmses))
    print("min_RMSE", min(rmses))
    print("median_RMSE", statistics.median(rmses))
    print("mean_RMSE", statistics.mean(rmses))
    print("rmse", statistics.median(rmses))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('--conformation_prediction')
    parser.add_argument('--conformation_answerkey')
    parser.add_argument('--molWeight_prediction')
    parser.add_argument('--molWeight_answerkey')
    parser.add_argument('--scores')
    

    args = parser.parse_args()

    if args.command == 'score-evaluation':
        score_evaluation(args.conformation_prediction, args.conformation_answerkey,
                     args.molWeight_prediction, args.molWeight_answerkey)

    elif args.command == 'score-submissionrun':
        score_submissionrun(args.scores)

    else:
        print("Invalid command: choose from \'score-evaluation\' or \'score-submissionrun\'.")
