"""
Wrapper to abstract away differences between Docker and Singularity
"""


from .engines import docker_engine
from .engines import singularity_engine
