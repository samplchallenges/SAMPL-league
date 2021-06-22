from setuptools import setup

setup(
    name='subtract',
    version='0.1',
    py_modules=['subtract'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        score=subtract:cli
    ''',
)
