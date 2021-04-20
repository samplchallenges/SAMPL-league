from setuptools import setup

setup(
    name='LogP-Calc',
    version='0.1',
    py_modules=['print_logP'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        print-LogP=print_logP:get_LogP
    ''',
)
