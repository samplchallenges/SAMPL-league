from .utils import GUEST_OUTPUT_DIR

from .docker_engine import DockerEngine
from .singularity_engine import SingularityEngine


REGISTERED_ENGINES = {}


def _register_engine(engine):
    REGISTERED_ENGINES[engine.name] = engine


_register_engine(DockerEngine)
_register_engine(SingularityEngine)
