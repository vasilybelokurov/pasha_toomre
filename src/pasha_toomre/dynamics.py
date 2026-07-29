"""Equations of motion, initial conditions, and scalar event functions."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

from .config import SimulationConfig
from .errors import DynamicsError
from .softening import get_softening

State = NDArray[np.float64]
RhsFunction = Callable[[float, State], State]


def initial_state(config: SimulationConfig) -> State:
    """Construct the documented initially circular, Sun-comoving planet state."""

    sun_velocity = np.sqrt(1.0 / (2.0 * abs(config.Z0)))
    planet_azimuthal_velocity = np.sqrt(1.0 / config.r0)
    return np.array(
        [
            config.Z0,
            sun_velocity,
            config.r0,
            0.0,
            config.Z0,
            0.0,
            planet_azimuthal_velocity,
            sun_velocity,
        ],
        dtype=float,
    )


def _validate_state(state: State) -> None:
    if np.shape(state) != (8,):
        raise DynamicsError(f"Expected an eight-component state, got {np.shape(state)}")
    if not np.all(np.isfinite(state)):
        raise DynamicsError("The dynamical state contains a non-finite value")


def make_rhs(config: SimulationConfig) -> RhsFunction:
    """Build the autonomous ODE right-hand side for one configuration."""

    softening = get_softening(config.softening)
    collision_distance2 = config.collision_tolerance**2

    def rhs(_t: float, state: State) -> State:
        _validate_state(state)
        Z, Vz, x, y, z, vx, vy, vz = state

        sun_distance2 = x**2 + y**2 + (z - Z) ** 2
        intruder_distance2 = x**2 + y**2 + (z + Z) ** 2
        if sun_distance2 <= collision_distance2:
            raise DynamicsError("The planet reached the Sun collision threshold")
        if intruder_distance2 <= collision_distance2:
            raise DynamicsError("The planet reached the intruder collision threshold")

        sun_distance3 = sun_distance2**1.5
        intruder_distance3 = intruder_distance2**1.5
        ax = -x / sun_distance3 - x / intruder_distance3
        ay = -y / sun_distance3 - y / intruder_distance3
        az = -(z - Z) / sun_distance3 - (z + Z) / intruder_distance3

        return np.array(
            [Vz, softening(Z, config.R), vx, vy, vz, ax, ay, az],
            dtype=float,
        )

    return rhs


def heliocentric_radial_velocity(state: State) -> float:
    """Return the planet's three-dimensional radial velocity from the Sun."""

    _validate_state(state)
    Z, Vz, x, y, z, vx, vy, vz = state
    dzeta = z - Z
    distance = np.sqrt(x**2 + y**2 + dzeta**2)
    if distance == 0.0:
        raise DynamicsError("Heliocentric radial velocity is singular at r=0")
    return float((x * vx + y * vy + dzeta * (vz - Vz)) / distance)


def sun_crosses_zero(_t: float, state: State) -> float:
    """Return the scalar Sun coordinate used to locate the upward crossing."""

    return float(state[0])


sun_crosses_zero.terminal = True
sun_crosses_zero.direction = 1.0


def find_pericenter(_t: float, state: State) -> float:
    """Locate negative-to-positive heliocentric radial-velocity crossings."""

    return heliocentric_radial_velocity(state)


find_pericenter.terminal = False
find_pericenter.direction = 1.0


def find_apocenter(_t: float, state: State) -> float:
    """Locate positive-to-negative heliocentric radial-velocity crossings."""

    return heliocentric_radial_velocity(state)


find_apocenter.terminal = False
find_apocenter.direction = -1.0
