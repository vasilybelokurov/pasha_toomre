"""Command-line entry point and complete simulation workflow."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

from .config import SimulationConfig, earth_config, venus_config
from .diagnostics import (
    ConservationCheck,
    Extremum,
    angular_momentum_conservation,
    compute_extrema,
    diagnostics_for_tau,
)
from .errors import PashaToomreError
from .plotting import generate_plots
from .reporting import write_reports, write_summary_json
from .simulation import SimulationResult, build_report_tau_grid, run_simulation


@dataclass(frozen=True)
class RunArtifacts:
    """Successful workflow products returned to programmatic callers."""

    result: SimulationResult
    extrema: tuple[Extremum, ...]
    conservation: ConservationCheck
    files: dict[str, Path]


def execute(config: SimulationConfig) -> RunArtifacts:
    """Run one simulation and generate all requested artifacts."""

    result = run_simulation(config)
    report_tau = build_report_tau_grid(result)
    diagnostics = diagnostics_for_tau(result, report_tau)
    extrema = compute_extrema(result)
    conservation = angular_momentum_conservation(result)
    files = write_reports(result, diagnostics, extrema, conservation)
    files.update(generate_plots(result))
    summary_path = config.output_dir / "run_summary.json"
    files["summary_json"] = summary_path
    write_summary_json(summary_path, result, extrema, conservation, files)
    return RunArtifacts(
        result=result,
        extrema=tuple(extrema),
        conservation=conservation,
        files=files,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a reconstructed Sun-intruder planetary encounter"
    )
    parser.add_argument("--planet", choices=("venus", "earth"), default="venus")
    parser.add_argument("--layout", choices=("overview", "extended"), default="overview")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--tau-final", type=float)
    parser.add_argument("--rtol", type=float, help="Override main relative tolerance")
    parser.add_argument("--atol", type=float, help="Override main absolute tolerance")
    parser.add_argument("--plot-dpi", type=int)
    parser.add_argument("--show", action="store_true", help="Request interactive display")
    return parser


def config_from_args(args: argparse.Namespace) -> SimulationConfig:
    """Build a validated configuration from parsed CLI arguments."""

    factory = venus_config if args.planet == "venus" else earth_config
    config = factory(layout=args.layout, output_dir=args.output_dir)
    overrides = {}
    if args.tau_final is not None:
        overrides["tau_final"] = args.tau_final
    if args.rtol is not None:
        overrides["main_rtol"] = args.rtol
    if args.atol is not None:
        overrides["main_atol"] = args.atol
    if args.plot_dpi is not None:
        overrides["plot_dpi"] = args.plot_dpi
    if args.show:
        overrides["show_plots"] = True
    return replace(config, **overrides)


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point. Return a conventional process exit status."""

    args = _parser().parse_args(argv)
    try:
        artifacts = execute(config_from_args(args))
    except (PashaToomreError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    result = artifacts.result
    check = artifacts.conservation
    print(f"Completed {result.config.run_name} in {result.config.output_dir}")
    print(f"Z=0 crossing: t_c={result.t_c:.12f}, tau_start={result.tau_start:.8f}")
    print(
        f"Lz: initial={check.initial_lz:.15f}, final={check.final_lz:.15f}, "
        f"|delta|={check.absolute_delta_lz:.3e}"
    )
    print(f"Detected {len(artifacts.extrema)} radial extrema")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
