from setuptools import setup

setup(
    name='OE-Docking',
    version='0.1',
    py_modules=['oedock_process'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        oedock=oedock_process:dock
    ''',
)
