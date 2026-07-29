"""Derived orbital quantities, event extrema, and conservation checks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .config import SimulationConfig
from .errors import DiagnosticError
from .simulation import SimulationResult

BASE_DIAGNOSTICS = (
    "Z",
    "zeta",
    "r",
    "rxy",
    "phi_phase",
    "vr",
    "vphi",
    "vz_rel",
    "vr3d",
    "e",
    "ra",
    "inclination",
)

EXTENDED_DIAGNOSTICS = (
    "Z",
    "zeta",
    "eta",
    "r",
    "rxy",
    "phi_phase",
    "vr",
    "vphi",
    "vz_rel",
    "vr3d",
    "e",
    "ra",
    "inclination",
)


@dataclass(frozen=True)
class Extremum:
    """One event-derived heliocentric radial extremum."""

    kind: str
    tau: float
    radius: float


@dataclass(frozen=True)
class ConservationCheck:
    """Initial and final axial angular momentum."""

    initial_lz: float
    final_lz: float
    absolute_delta_lz: float


def diagnostic_names(layout: str) -> tuple[str, ...]:
    """Return the ordered report/plot diagnostics for a layout."""

    if layout == "overview":
        return BASE_DIAGNOSTICS
    if layout == "extended":
        return EXTENDED_DIAGNOSTICS
    raise DiagnosticError(f"Unknown diagnostic layout: {layout!r}")


def compute_diagnostics(
    states: ArrayLike,
    tau: ArrayLike,
    config: SimulationConfig,
) -> dict[str, NDArray[np.float64]]:
    """Compute all common and extended diagnostics for vectorized states."""

    state_values = np.asarray(states, dtype=float)
    tau_values = np.atleast_1d(np.asarray(tau, dtype=float))
    if state_values.ndim == 1:
        state_values = state_values[:, np.newaxis]
    if state_values.shape != (8, tau_values.size):
        raise DiagnosticError(
            f"Expected states with shape (8, {tau_values.size}), got {state_values.shape}"
        )
    if not np.all(np.isfinite(state_values)) or not np.all(np.isfinite(tau_values)):
        raise DiagnosticError("States and tau values must be finite")

    Z, Vz, x, y, z, vx, vy, vz = state_values
    zeta = z - Z
    eta = z + Z
    rxy = np.sqrt(x**2 + y**2)
    radius = np.sqrt(rxy**2 + zeta**2)
    if np.any(radius <= config.collision_tolerance):
        raise DiagnosticError("A heliocentric diagnostic reached the collision threshold")

    phi_phase = np.mod(np.arctan2(y, x) / (2.0 * np.pi), 1.0)
    vr = np.zeros_like(radius)
    vphi = np.zeros_like(radius)
    planar = rxy > config.collision_tolerance
    vr[planar] = (x[planar] * vx[planar] + y[planar] * vy[planar]) / rxy[planar]
    vphi[planar] = (
        x[planar] * vy[planar] - y[planar] * vx[planar]
    ) / rxy[planar]

    vz_rel = vz - Vz
    vr3d = (x * vx + y * vy + zeta * vz_rel) / radius
    Lx = y * vz_rel - zeta * vy
    Ly = zeta * vx - x * vz_rel
    Lz = x * vy - y * vx
    L2 = Lx**2 + Ly**2 + Lz**2
    L = np.sqrt(L2)
    inclination = np.degrees(np.arctan2(np.sqrt(Lx**2 + Ly**2), Lz))

    eccentricity = np.full_like(radius, np.nan)
    aphelion = np.full_like(radius, np.nan)
    active = tau_values >= config.eccentricity_start_tau
    if np.any(active):
        speed2 = vx[active] ** 2 + vy[active] ** 2 + vz_rel[active] ** 2
        radicand = 1.0 + (speed2 - 2.0 / radius[active]) * L2[active]
        minimum = float(np.min(radicand))
        if minimum < -config.eccentricity_tolerance:
            raise DiagnosticError(
                "Osculating eccentricity has a materially negative radicand: "
                f"{minimum:.6e}"
            )
        radicand[np.abs(radicand) < config.eccentricity_tolerance] = 0.0
        radicand = np.maximum(radicand, 0.0)
        eccentricity[active] = np.sqrt(radicand)
        bound = active & (eccentricity < 1.0)
        aphelion[bound] = L2[bound] / (1.0 - eccentricity[bound])

    values = {
        "tau": tau_values,
        "Z": Z,
        "zeta": zeta,
        "eta": eta,
        "r": radius,
        "rxy": rxy,
        "phi_phase": phi_phase,
        "vr": vr,
        "vphi": vphi,
        "vz_rel": vz_rel,
        "vr3d": vr3d,
        "e": eccentricity,
        "ra": aphelion,
        "inclination": inclination,
        "Lx": Lx,
        "Ly": Ly,
        "Lz": Lz,
        "L": L,
    }
    for name, array in values.items():
        if array.shape != tau_values.shape:
            raise DiagnosticError(f"Diagnostic {name!r} has an inconsistent shape")
        if name not in ("e", "ra") and not np.all(np.isfinite(array)):
            raise DiagnosticError(f"Diagnostic {name!r} contains a non-finite value")
    if np.any(~np.isfinite(eccentricity[active])):
        raise DiagnosticError("Active eccentricity values must be finite")
    return values


def diagnostics_for_tau(
    result: SimulationResult, tau: ArrayLike
) -> dict[str, NDArray[np.float64]]:
    """Evaluate the dense solution and compute diagnostics on a tau grid."""

    tau_values = np.atleast_1d(np.asarray(tau, dtype=float))
    return compute_diagnostics(result.evaluate_tau(tau_values), tau_values, result.config)


def compute_extrema(result: SimulationResult) -> list[Extremum]:
    """Convert event roots into sorted heliocentric radial extrema."""

    extrema: list[Extremum] = []
    event_groups = (
        ("pericenter", result.solution.t_events[0]),
        ("apocenter", result.solution.t_events[1]),
    )
    for kind, event_times in event_groups:
        for event_time in event_times:
            state = np.asarray(result.solution.sol(float(event_time)), dtype=float)
            Z, _, x, y, z, _, _, _ = state
            radius = float(np.sqrt(x**2 + y**2 + (z - Z) ** 2))
            tau = float(result.time_to_tau(float(event_time)))
            extrema.append(Extremum(kind=kind, tau=tau, radius=radius))
    extrema.sort(key=lambda item: item.tau)
    return extrema


def angular_momentum_conservation(result: SimulationResult) -> ConservationCheck:
    """Return the initial/final heliocentric axial angular-momentum check."""

    start = np.asarray(result.solution.sol(0.0), dtype=float)
    final = np.asarray(result.solution.sol(result.t_end), dtype=float)
    initial_lz = float(start[2] * start[6] - start[3] * start[5])
    final_lz = float(final[2] * final[6] - final[3] * final[5])
    return ConservationCheck(
        initial_lz=initial_lz,
        final_lz=final_lz,
        absolute_delta_lz=abs(final_lz - initial_lz),
    )
