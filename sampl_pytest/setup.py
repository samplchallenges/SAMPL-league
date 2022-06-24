from setuptools import setup

setup(
    name="sampl-pytest",
    packages=["sampl_pytest"],
    # the following makes a plugin available to pytest
    entry_points={"pytest11": ["name_of_plugin = sampl_pytest.hook_config"]},
    # custom PyPI classifier for pytest plugins
    classifiers=["Framework :: Pytest"],
)