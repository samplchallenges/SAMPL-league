import math

import click

# --predicted_molwt = 78
# --expected_molwt = 78.2156
# --predicted_
#@click.command()
#@click.argument("a", type=click.FLOAT)
#@click.argument("b", type=click.FLOAT)
#def subtract_floats(a, b):
#    print(a-b)

@click.group()
def cli():
    pass

@cli.command()
@click.option("--item", nargs=3, type=(str, click.FLOAT, click.FLOAT), multiple=True)
def score_evaluation(item):
    # ignore error_estimate
    for l_item in item:
        key, expected_value, predicted_value = l_item
        if key == 'molwt':
            print (expected_value - predicted_value)


@cli.command()
@click.argument("scores", nargs=-1, type=click.FLOAT)
def score_submissionrun(scores):
    print("RMSE", math.sqrt(sum(score*score for score in scores)/len(scores)))
    print("MAU", sum(abs(score) for score in scores)/len(scores))
