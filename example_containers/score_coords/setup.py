from setuptools import setup

setup(
    name='score_coords',
    version='0.1',
    py_modules=['score'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        score-coords=score:score_coords
    ''',
)
