from setuptools import setup

setup(
    name='logging_example',
    version='0.1',
    py_modules=['logging_example'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        calc-coords=logging_example:calc_coords
    ''',
)
