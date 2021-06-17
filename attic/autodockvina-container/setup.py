from setuptools import setup

setup(
    name='Autodock',
    version='0.1',
    py_modules=['run_autodock'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        RunAutodock=run_autodock:autodock
    ''',
)
