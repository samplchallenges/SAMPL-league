from setuptools import setup

setup(
    name='AutoDock-rdkit',
    version='0.1',
    py_modules=[
		'run_autodock',
	],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        run-autodock=run_autodock:autodock
    ''',
)
