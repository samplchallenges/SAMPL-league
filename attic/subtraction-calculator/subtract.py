import click


@click.command()
@click.argument("a", type=click.FLOAT)
@click.argument("b", type=click.FLOAT)
def subtract_floats(a, b):
    print(a-b)
