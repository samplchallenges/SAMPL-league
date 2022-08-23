from email import header
import json
import os
import os.path
import random
import statistics
import argparse
import csv

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

def _load_csv(filename):
    result_dict = {}
    with open(filename, "r", encoding="utf8") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            result_dict[row["id"]] = float(row["value"])
    return result_dict


def _write_csv(filename, score_dict):
    with open(filename, "w", encoding="utf8") as fp:
        writer = csv.DictWriter(fp, ["id", "value"])
        writer.writeheader()
        for k, v in score_dict.items():
            writer.writerow({"id": k, "value": v})


def score_batch(output_dir, conformation_prediction, conformation_answerkey,
                molweight_prediction, molweight_answerkey):
    # Files
    # ignore conformation for now
    if output_dir is None:
        output_dir = "/mnt/outputs"
    score_dict = {}
    molweight_predictions = _load_csv(molweight_prediction)
    molweight_answers = _load_csv(molweight_answerkey)
    if len(molweight_answers) != len(molweight_predictions):
        raise ValueError("Predictions and Answers don't match {len(molweight_predictions)} vs {len(molweight_answers)}")
    for key, prediction in molweight_predictions.items():
        try:
            answer = molweight_answers[key]
            score_dict[key] = abs(prediction - answer)
        except KeyError:
            raise ValueError("No answer found for key {key}")
    rmse_path = os.path.join(output_dir, "rmse.csv")
    diff_path = os.path.join(output_dir, "diff.csv")
    _write_csv(rmse_path, score_dict)
    _write_csv(diff_path, score_dict)
    print("rmse", rmse_path)
    print("diff", diff_path)


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
    parser.add_argument('--output-dir', help="Output Dir (for batching)")
 

    args = parser.parse_args()

    if args.command == 'score-evaluation':
        score_evaluation(args.conformation_prediction, args.conformation_answerkey,
                         args.molWeight_prediction, args.molWeight_answerkey)
    elif args.command == 'score-batch':
        score_batch(args.output_dir, args.conformation_prediction, args.conformation_answerkey,
                    args.molWeight_prediction, args.molWeight_answerkey)
    elif args.command == 'score-submissionrun':
        score_submissionrun(args.scores)
    else:
        print("Invalid command: choose from \'score-evaluation\', \'score-batch\', or \'score-submissionrun\'.")
