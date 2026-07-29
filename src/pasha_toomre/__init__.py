"""Sun-intruder encounter simulation library."""

from .config import SimulationConfig, earth_config, venus_config
from .simulation import SimulationResult, run_simulation

__all__ = [
    "SimulationConfig",
    "SimulationResult",
    "earth_config",
    "run_simulation",
    "venus_config",
]

__version__ = "0.1.0"
