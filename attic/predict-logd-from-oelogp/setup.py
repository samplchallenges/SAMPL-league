from setuptools import setup

setup(
    name='LogD-Calc',
    version='0.0',
    py_modules=['print_logd'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        print-LogD=print_logd:get_logd
    ''',
)
