"""Two-pass integration and time-grid generation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from .config import SimulationConfig
from .dynamics import (
    find_apocenter,
    find_pericenter,
    initial_state,
    make_rhs,
    sun_crosses_zero,
)
from .errors import IntegrationError


@dataclass(frozen=True)
class SimulationResult:
    """Successful probe and main integrations plus their time transformation."""

    config: SimulationConfig
    t_c: float
    tau_start: float
    t_end: float
    probe: Any
    solution: Any

    def tau_to_time(self, tau: ArrayLike) -> NDArray[np.float64]:
        values = np.asarray(tau, dtype=float)
        return values * np.pi / 6.0 + self.t_c

    def time_to_tau(self, time: ArrayLike) -> NDArray[np.float64]:
        values = np.asarray(time, dtype=float)
        return 6.0 * (values - self.t_c) / np.pi

    def evaluate_tau(self, tau: ArrayLike) -> NDArray[np.float64]:
        values = np.atleast_1d(np.asarray(tau, dtype=float))
        tolerance = 1.0e-10
        if np.any(values < self.tau_start - tolerance):
            raise IntegrationError(
                f"Requested tau below integrated start {self.tau_start:.8f}"
            )
        if np.any(values > self.config.tau_final + tolerance):
            raise IntegrationError(
                f"Requested tau above integrated end {self.config.tau_final:.8f}"
            )
        states = np.asarray(self.solution.sol(self.tau_to_time(values)), dtype=float)
        if states.shape != (8, values.size) or not np.all(np.isfinite(states)):
            raise IntegrationError("Dense output returned invalid state data")
        return states


def run_simulation(config: SimulationConfig) -> SimulationResult:
    """Run the crossing probe and full dense-output encounter integration."""

    rhs = make_rhs(config)
    state0 = initial_state(config)
    probe_end = config.tau_final * np.pi / 6.0 + 500.0
    probe = solve_ivp(
        fun=rhs,
        t_span=(0.0, probe_end),
        y0=state0,
        method=config.method,
        rtol=config.probe_rtol,
        atol=config.probe_atol,
        events=sun_crosses_zero,
    )
    if not probe.success:
        raise IntegrationError(f"Probe integration failed: {probe.message}")
    if len(probe.t_events) != 1 or probe.t_events[0].size != 1:
        raise IntegrationError("The probe did not find exactly one upward Z=0 crossing")
    t_c = float(probe.t_events[0][0])
    crossing_state = np.asarray(probe.y_events[0][0], dtype=float)
    if abs(crossing_state[0]) > 1.0e-8:
        raise IntegrationError(
            f"Detected crossing does not satisfy Z=0: Z={crossing_state[0]:.3e}"
        )

    t_end = t_c + config.tau_final * np.pi / 6.0
    solution = solve_ivp(
        fun=rhs,
        t_span=(0.0, t_end),
        y0=state0,
        method=config.method,
        rtol=config.main_rtol,
        atol=config.main_atol,
        events=(find_pericenter, find_apocenter),
        dense_output=True,
    )
    if not solution.success:
        raise IntegrationError(f"Main integration failed: {solution.message}")
    if solution.sol is None:
        raise IntegrationError("Main integration did not return dense output")
    if not np.all(np.isfinite(solution.y)):
        raise IntegrationError("Main integration contains non-finite state values")

    tau_start = -6.0 * t_c / np.pi
    result = SimulationResult(
        config=config,
        t_c=t_c,
        tau_start=tau_start,
        t_end=t_end,
        probe=probe,
        solution=solution,
    )
    result.evaluate_tau(np.array([tau_start, config.tau_final]))
    return result


def build_report_tau_grid(result: SimulationResult) -> NDArray[np.float64]:
    """Reproduce the progressively refined monthly grid from the documents."""

    tau_current = float(np.ceil(result.tau_start))
    nodes: list[float] = []
    segments = [
        (-3.0, 1.0),
        (-1.0, 0.1),
        (-0.5, 0.05),
        (-0.2, 0.02),
        (-0.1, 0.01),
        (-0.01, 0.005),
        (-0.002, 0.002),
        (0.002, 0.0004),
        (0.01, 0.002),
        (0.1, 0.005),
        (0.2, 0.01),
        (0.5, 0.02),
        (1.0, 0.05),
        (4.0, 0.1),
        (result.config.tau_final, 0.5),
    ]
    for endpoint, step in segments:
        upper = min(endpoint, result.config.tau_final)
        while tau_current <= upper + 1.0e-12:
            if tau_current >= result.tau_start - 1.0e-12:
                nodes.append(tau_current)
            tau_current += step
        if endpoint >= result.config.tau_final:
            break

    rounded = {
        round(value, 5)
        for value in nodes
        if result.tau_start <= value <= result.config.tau_final
    }
    rounded.add(0.0)
    rounded.add(round(result.config.tau_final, 5))
    if not rounded:
        raise IntegrationError("The adaptive report grid is empty")
    return np.asarray(sorted(rounded), dtype=float)


def build_plot_tau_grid(
    result: SimulationResult, start: float, end: float, points: int
) -> NDArray[np.float64]:
    """Return a validated uniform plot grid inside the integrated interval."""

    if points < 2 or start >= end:
        raise IntegrationError("Invalid plot-grid request")
    if start < result.tau_start - 1.0e-10:
        raise IntegrationError(
            f"Plot start {start} precedes integrated tau {result.tau_start:.8f}"
        )
    if end > result.config.tau_final + 1.0e-10:
        raise IntegrationError(
            f"Plot end {end} exceeds integrated tau {result.config.tau_final:.8f}"
        )
    return np.linspace(start, end, points)
