from setuptools import setup

setup(
    name='dock',
    version='0.1',
    py_modules=['run_autodock'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        dock=run_autodock:autodock
    ''',
)
