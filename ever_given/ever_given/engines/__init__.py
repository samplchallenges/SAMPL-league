from .utils import GUEST_OUTPUT_DIR

from .docker_engine import DockerEngine
from .singularity_engine import SingularityEngine



REGISTERED_ENGINES = {}


def _register_engine(engine_class):
    REGISTERED_ENGINES[engine_class.name()] = engine_class


_register_engine(DockerEngine)
_register_engine(SingularityEngine)
