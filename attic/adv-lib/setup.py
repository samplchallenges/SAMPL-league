from setuptools import setup

setup(
    name='Autodock-Vina',
    version='0.1',
    py_modules=['autodockvina'],
    install_requires=[
        'Click',
	'Vina',
    ],
    entry_points='''
        [console_scripts]
        adv=autodockvina:run_autodock
    ''',
)
