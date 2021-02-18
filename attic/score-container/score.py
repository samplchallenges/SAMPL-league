import click


@click.command()
@click.argument("path", type=click.Path(exists=True))
def score(path):
    from pathlib import Path

    p = Path(path)

    # Even if the sorting is not "right" it will be equally wrong
    predictions = sorted(list(p.glob("prediction_mol_*")))
    measurements = sorted(list(p.glob("measurement_mol_*")))

    assert len(predictions) == len(measurements)

    scores = []
    for prediction, measurement in zip(predictions, measurements):
        with open(prediction, "r") as prediction_file, open(
            measurement, "r"
        ) as measurement_file:
            prediction_val, measurement_val = (
                float(prediction_file.readline().strip()),
                float(measurement_file.readline().strip()),
            )
        score = (prediction_val - measurement_val) ** 2
        scores.append(score)
        mol_name = prediction.name[len("prediction_"):]
        score_file_name = prediction.parent.joinpath(f"score_{mol_name}")
        with open(score_file_name, "w") as score_file:
            score_file.write(f"{str(score)}\n")

    print(f"total RMSE: {sum(scores)}")


if __name__ == "__main__":
    score()
