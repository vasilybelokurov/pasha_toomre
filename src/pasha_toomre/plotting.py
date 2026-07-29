"""Reconstructed overview and extended plot layouts."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

_cache_root = Path(tempfile.gettempdir()) / "pasha_toomre_cache"
_mpl_cache = _cache_root / "matplotlib"
_xdg_cache = _cache_root / "xdg"
_mpl_cache.mkdir(parents=True, exist_ok=True)
_xdg_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_cache))
os.environ.setdefault("XDG_CACHE_HOME", str(_xdg_cache))

import matplotlib

matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
import matplotlib.pyplot as plt
import numpy as np

from .diagnostics import diagnostics_for_tau
from .simulation import SimulationResult, build_plot_tau_grid

LABELS = {
    "Z": r"$Z$",
    "zeta": r"$\zeta$",
    "eta": r"$\eta$",
    "r": r"$r$",
    "rxy": r"$r_{xy}$",
    "phi_phase": r"$\phi_{\rm phase}$",
    "vr": r"$v_r$",
    "vphi": r"$v_\phi$",
    "vz_rel": r"$v_{z,\rm rel}$",
    "vr3d": r"$v_{r,3D}$",
    "e": r"$e$",
    "ra": r"$r_a$",
    "inclination": r"$i\ (\rm deg)$",
}

OVERVIEW_NAMES = (
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

COORDINATE_NAMES = ("Z", "zeta", "eta", "r", "rxy", "phi_phase")
ORBIT_NAMES = ("vr", "vphi", "vz_rel", "vr3d", "e", "ra", "inclination")


def _draw_grid(
    tau: np.ndarray,
    diagnostics: dict[str, np.ndarray],
    names: tuple[str, ...],
    rows: int,
    columns: int,
    figure_size: tuple[float, float],
    title: str,
    color: str,
    output_path: Path,
    dpi: int,
    keep_open: bool = False,
) -> Path:
    figure, axes = plt.subplots(rows, columns, figsize=figure_size, sharex=True)
    flat_axes = np.atleast_1d(axes).ravel()
    for index, name in enumerate(names):
        axis = flat_axes[index]
        axis.plot(tau, diagnostics[name], color=color, linewidth=1.5, label="Trajectory")
        axis.axvline(0.0, color="red", linestyle=":", alpha=0.7, label=r"$\tau=0$")
        axis.set_ylabel(LABELS[name], fontsize=11)
        axis.grid(True, linestyle="--", alpha=0.45)
    for axis in flat_axes[len(names) :]:
        axis.set_visible(False)
    for column in range(columns):
        visible_indices = [
            index for index in range(len(names)) if index % columns == column
        ]
        if visible_indices:
            flat_axes[visible_indices[-1]].set_xlabel(r"$\tau$", fontsize=11)
    flat_axes[0].legend(loc="best", fontsize=8)
    figure.suptitle(title, fontsize=14, fontweight="bold", y=0.995)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.985))
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if output_path.stat().st_size == 0:
        raise RuntimeError(f"Generated an empty plot: {output_path}")
    if not keep_open:
        plt.close(figure)
    return output_path


def generate_plots(result: SimulationResult) -> dict[str, Path]:
    """Generate the selected document-variant plot layout."""

    config = result.config
    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    full_tau = build_plot_tau_grid(
        result, config.full_tau_start, config.tau_final, config.full_plot_points
    )
    detailed_tau = build_plot_tau_grid(
        result,
        config.detailed_tau_start,
        min(config.detailed_tau_end, config.tau_final),
        config.detailed_plot_points,
    )
    full_data = diagnostics_for_tau(result, full_tau)
    detailed_data = diagnostics_for_tau(result, detailed_tau)

    outputs: dict[str, Path] = {}
    if config.layout == "overview":
        outputs["plot_full_overview"] = _draw_grid(
            full_tau,
            full_data,
            OVERVIEW_NAMES,
            6,
            2,
            (14.0, 20.0),
            f"Trajectory Parameters: {config.planet_name} "
            rf"($\tau\in[{config.full_tau_start:g},{config.tau_final:g}]$)",
            "blue",
            output_dir / "trajectory_full_overview.png",
            config.plot_dpi,
            config.show_plots,
        )
        outputs["plot_singularity_core"] = _draw_grid(
            detailed_tau,
            detailed_data,
            OVERVIEW_NAMES,
            6,
            2,
            (14.0, 20.0),
            f"Detailed Encounter: {config.planet_name} "
            rf"($\tau\in[{config.detailed_tau_start:g},{min(config.detailed_tau_end, config.tau_final):g}]$)",
            "darkgreen",
            output_dir / "trajectory_singularity_core.png",
            config.plot_dpi,
            config.show_plots,
        )
    else:
        plot_specs = (
            (
                "plot_full_coordinates",
                full_tau,
                full_data,
                COORDINATE_NAMES,
                3,
                (14.0, 10.0),
                "Full Range: Coordinates and Angles",
                "blue",
                "trajectory_full_coordinates.png",
            ),
            (
                "plot_full_orbits",
                full_tau,
                full_data,
                ORBIT_NAMES,
                4,
                (14.0, 13.0),
                "Full Range: Velocities and Orbital Parameters",
                "blue",
                "trajectory_full_velocities.png",
            ),
            (
                "plot_detailed_coordinates",
                detailed_tau,
                detailed_data,
                COORDINATE_NAMES,
                3,
                (14.0, 10.0),
                "Central Encounter: Coordinates and Angles",
                "darkgreen",
                "trajectory_detailed_coordinates.png",
            ),
            (
                "plot_detailed_orbits",
                detailed_tau,
                detailed_data,
                ORBIT_NAMES,
                4,
                (14.0, 13.0),
                "Central Encounter: Velocities and Orbital Parameters",
                "darkgreen",
                "trajectory_detailed_velocities.png",
            ),
        )
        for key, tau, data, names, rows, size, title, color, filename in plot_specs:
            outputs[key] = _draw_grid(
                tau,
                data,
                names,
                rows,
                2,
                size,
                f"{title}: {config.planet_name}",
                color,
                output_dir / filename,
                config.plot_dpi,
                config.show_plots,
            )
    if config.show_plots:
        plt.show()
        plt.close("all")
    return outputs
