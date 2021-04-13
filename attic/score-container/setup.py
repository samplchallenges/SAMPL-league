from setuptools import setup

setup(
    name='score-submission',
    version='0.1',
    py_modules=['score'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        score-submission=score:score
    ''',
)
