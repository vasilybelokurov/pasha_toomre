"""Configuration objects and documented planet presets."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from .errors import ConfigurationError

Layout = Literal["overview", "extended"]


@dataclass(frozen=True)
class SimulationConfig:
    """Complete, immutable configuration for one encounter simulation."""

    run_name: str
    planet_name: str
    layout: Layout
    output_dir: Path
    r0: float = 0.7226
    Z0: float = -5.0
    R: float = 0.00464913034
    tau_final: float = 80.0
    softening: str = "plummer4"
    method: str = "DOP853"
    probe_rtol: float = 1.0e-11
    probe_atol: float = 1.0e-13
    main_rtol: float = 1.0e-12
    main_atol: float = 1.0e-14
    full_tau_start: float = -10.0
    detailed_tau_start: float = -2.0
    detailed_tau_end: float = 4.0
    eccentricity_start_tau: float = 0.1
    eccentricity_tolerance: float = 1.0e-10
    collision_tolerance: float = 1.0e-12
    full_plot_points: int = 2000
    detailed_plot_points: int = 1000
    plot_dpi: int = 200
    show_plots: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "output_dir", Path(self.output_dir))
        if self.layout not in ("overview", "extended"):
            raise ConfigurationError(f"Unknown layout: {self.layout!r}")
        if not self.run_name.strip() or not self.planet_name.strip():
            raise ConfigurationError("run_name and planet_name must be non-empty")
        if self.r0 <= 0.0:
            raise ConfigurationError("r0 must be positive")
        if self.Z0 >= 0.0:
            raise ConfigurationError("Z0 must be negative for an upward crossing")
        if self.R <= 0.0:
            raise ConfigurationError("R must be positive")
        if self.tau_final <= 0.0:
            raise ConfigurationError("tau_final must be positive")
        if self.detailed_tau_start >= self.detailed_tau_end:
            raise ConfigurationError("detailed plot bounds are reversed")
        if self.full_tau_start >= self.tau_final:
            raise ConfigurationError("full plot start must precede tau_final")
        for name in ("probe_rtol", "probe_atol", "main_rtol", "main_atol"):
            if getattr(self, name) <= 0.0:
                raise ConfigurationError(f"{name} must be positive")
        if self.eccentricity_tolerance <= 0.0 or self.collision_tolerance <= 0.0:
            raise ConfigurationError("diagnostic tolerances must be positive")
        if self.full_plot_points < 2 or self.detailed_plot_points < 2:
            raise ConfigurationError("plot grids need at least two points")
        if self.plot_dpi <= 0:
            raise ConfigurationError("plot_dpi must be positive")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe configuration mapping."""

        values = asdict(self)
        values["output_dir"] = str(self.output_dir)
        return values


def venus_config(
    layout: Layout = "overview", output_dir: str | Path | None = None
) -> SimulationConfig:
    """Return the documented Venus configuration."""

    run_name = "test_venus" if layout == "overview" else "venus_extended"
    destination = Path(output_dir) if output_dir is not None else Path("outputs") / run_name
    return SimulationConfig(
        run_name=run_name,
        planet_name="Venus",
        layout=layout,
        output_dir=destination,
        r0=0.7226,
    )


def earth_config(
    layout: Layout = "overview", output_dir: str | Path | None = None
) -> SimulationConfig:
    """Return an Earth-radius configuration for later reference studies."""

    run_name = "earth_overview" if layout == "overview" else "earth_extended"
    destination = Path(output_dir) if output_dir is not None else Path("outputs") / run_name
    return SimulationConfig(
        run_name=run_name,
        planet_name="Earth",
        layout=layout,
        output_dir=destination,
        r0=1.0,
    )
