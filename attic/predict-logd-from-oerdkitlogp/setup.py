from setuptools import setup

setup(
    name='oechem-rdkit-logd',
    version='0.1',
    py_modules=['print_logD'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        print-LogD=print_logD:get_logd
    ''',
)
