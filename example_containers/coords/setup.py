from setuptools import setup

setup(
    name='coords',
    version='0.1',
    py_modules=['coords'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        calc-coords=coords:calc_coords
    ''',
)
