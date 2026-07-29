"""Interchangeable Sun-intruder gravitational-softening laws."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .errors import ConfigurationError, UnknownSofteningError

SofteningFunction = Callable[[ArrayLike, float], float | NDArray[np.float64]]


def plummer4(Z: ArrayLike, R: float) -> float | NDArray[np.float64]:
    """Return the fourth-order Plummer acceleration documented in the files."""

    if R <= 0.0:
        raise ConfigurationError("R must be positive")
    coordinates = np.asarray(Z, dtype=float)
    q = coordinates / R
    acceleration = -(
        q * (q**2 + 5.0 / 8.0)
    ) / (4.0 * R**2 * (q**2 + 1.0 / 4.0) ** 2.5)
    if acceleration.ndim == 0:
        return float(acceleration)
    return acceleration


SOFTENING_LAWS: dict[str, SofteningFunction] = {
    "plummer4": plummer4,
}


def get_softening(name: str) -> SofteningFunction:
    """Resolve a registered softening law by name."""

    try:
        return SOFTENING_LAWS[name]
    except KeyError as exc:
        available = ", ".join(sorted(SOFTENING_LAWS))
        raise UnknownSofteningError(
            f"Unknown softening law {name!r}; available laws: {available}"
        ) from exc
