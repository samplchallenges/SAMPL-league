
# Allow command-line test runs to hit just docker, just singularity, or both
def pytest_addoption(parser):
    parser.addoption("--docker", dest="engines", action="append_const", const="docker")
    parser.addoption(
        "--singularity", dest="engines", action="append_const", const="singularity"
    )


def pytest_generate_tests(metafunc):
    if "container_engine" in metafunc.fixturenames:
        engines = metafunc.config.getoption("engines")
        metafunc.parametrize("container_engine", engines)
