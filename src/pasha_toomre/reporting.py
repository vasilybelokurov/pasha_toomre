"""Human-readable and machine-readable simulation outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .diagnostics import ConservationCheck, Extremum, diagnostic_names
from .simulation import SimulationResult


def _table_columns(layout: str) -> tuple[str, ...]:
    return ("tau",) + diagnostic_names(layout)


def write_trajectory_csv(
    destination: Path,
    diagnostics: dict[str, np.ndarray],
    layout: str,
) -> Path:
    """Write the diagnostic table as numeric CSV."""

    columns = _table_columns(layout)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for index in range(diagnostics["tau"].size):
            writer.writerow([f"{float(diagnostics[name][index]):.12g}" for name in columns])
    return destination


def _display_value(name: str, value: float) -> str:
    if not np.isfinite(value):
        return "---" if name in ("e", "ra") else "nan"
    if name == "tau":
        return f"{value:.5f}"
    if name == "inclination":
        return f"{value:.4f}"
    return f"{value:.8f}"


def write_trajectory_text(
    destination: Path,
    result: SimulationResult,
    diagnostics: dict[str, np.ndarray],
    layout: str,
    extrema: list[Extremum],
    conservation: ConservationCheck,
) -> Path:
    """Write a readable reconstruction of the original fixed-width report."""

    columns = _table_columns(layout)
    widths = {name: max(12, len(name) + 2) for name in columns}
    with destination.open("w", encoding="utf-8") as handle:
        handle.write(
            f"ORBITAL SIMULATION REPORT: {result.config.planet_name.upper()}\n"
            f"run={result.config.run_name}, layout={layout}, "
            f"softening={result.config.softening}\n"
            f"r0={result.config.r0:.8f}, Z0={result.config.Z0:.8f}, "
            f"R={result.config.R:.11f}\n"
            f"Z=0 crossing t_c={result.t_c:.12f}; "
            f"tau_start={result.tau_start:.8f}; tau_final={result.config.tau_final:.4f}\n\n"
        )
        header = " | ".join(name.ljust(widths[name]) for name in columns)
        separator = "-+-".join("-" * widths[name] for name in columns)
        handle.write(header + "\n")
        handle.write(separator + "\n")
        for index in range(diagnostics["tau"].size):
            row = " | ".join(
                _display_value(name, float(diagnostics[name][index])).ljust(widths[name])
                for name in columns
            )
            handle.write(row + "\n")

        handle.write("\nEVENT-DERIVED RADIAL EXTREMA\n")
        handle.write("----------------------------\n")
        for extremum in extrema:
            handle.write(
                f"{extremum.kind:<10} tau={extremum.tau:.8f} "
                f"r={extremum.radius:.10f}\n"
            )
        handle.write("\nAXIAL ANGULAR-MOMENTUM CONSERVATION\n")
        handle.write("-----------------------------------\n")
        handle.write(f"initial Lz = {conservation.initial_lz:.15f}\n")
        handle.write(f"final Lz   = {conservation.final_lz:.15f}\n")
        handle.write(f"|delta Lz| = {conservation.absolute_delta_lz:.15e}\n")
    return destination


def write_extrema_csv(destination: Path, extrema: list[Extremum]) -> Path:
    """Write event-derived extrema as CSV."""

    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("kind", "tau", "radius"))
        for item in extrema:
            writer.writerow((item.kind, f"{item.tau:.12g}", f"{item.radius:.12g}"))
    return destination


def write_summary_json(
    destination: Path,
    result: SimulationResult,
    extrema: list[Extremum],
    conservation: ConservationCheck,
    output_files: dict[str, Path],
) -> Path:
    """Write configuration, solver, event, and output metadata."""

    payload: dict[str, Any] = {
        "config": result.config.to_dict(),
        "crossing": {
            "t_c": result.t_c,
            "tau_start": result.tau_start,
            "crossing_z": float(result.probe.y_events[0][0][0]),
        },
        "solver": {
            "method": result.config.method,
            "probe_success": bool(result.probe.success),
            "probe_nfev": int(result.probe.nfev),
            "main_success": bool(result.solution.success),
            "main_nfev": int(result.solution.nfev),
            "main_steps": int(result.solution.t.size),
        },
        "conservation": {
            "initial_lz": conservation.initial_lz,
            "final_lz": conservation.final_lz,
            "absolute_delta_lz": conservation.absolute_delta_lz,
        },
        "extrema": [
            {"kind": item.kind, "tau": item.tau, "radius": item.radius}
            for item in extrema
        ],
        "outputs": {name: str(path) for name, path in sorted(output_files.items())},
    }
    with destination.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return destination


def write_reports(
    result: SimulationResult,
    diagnostics: dict[str, np.ndarray],
    extrema: list[Extremum],
    conservation: ConservationCheck,
) -> dict[str, Path]:
    """Write all report files except the final run summary."""

    output_dir = result.config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "trajectory_text": write_trajectory_text(
            output_dir / "trajectory_report.txt",
            result,
            diagnostics,
            result.config.layout,
            extrema,
            conservation,
        ),
        "trajectory_csv": write_trajectory_csv(
            output_dir / "trajectory_report.csv", diagnostics, result.config.layout
        ),
        "extrema_csv": write_extrema_csv(output_dir / "extrema.csv", extrema),
    }
    return outputs
