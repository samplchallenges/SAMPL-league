from setuptools import setup

setup(
    name='molwt',
    version='0.1',
    py_modules=['molwt'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        calc-mw=molwt:calc_mol_wt
    ''',
)
