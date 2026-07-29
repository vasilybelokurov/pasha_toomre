# Pasha–Toomre Sun–Intruder Simulations

This repository reconstructs the two Venus programs supplied as Word documents. A shared Python library now contains the dynamics, integration, diagnostics, reporting, and plotting logic. Two thin drivers reproduce the original 12-diagnostic layout and the extended 13-diagnostic layout.

The original files in `docs/` remain unchanged.

## Requirements

- Python 3.10 or newer
- NumPy
- SciPy
- Matplotlib

Create a project-local virtual environment and install the package from the project
root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

The commands below assume this environment is active. If your system reserves
`python` for Python 2, use `python3` instead.

The implementation was verified with:

- Python 3.13
- NumPy 1.26.4
- SciPy 1.16.1
- Matplotlib 3.10.3

## Run the documented Venus variants

From the project root:

```bash
python scripts/run_test_venus.py
python scripts/run_venus_extended.py
```

The first command produces the two 6-by-2 figures described by `Test Venus`. The second adds the intruder-relative coordinate `eta = z + Z` to the report and divides 13 diagnostics across four figures.

The general interface is also available through:

```bash
python -m pasha_toomre.cli --planet venus --layout overview
python -m pasha_toomre.cli --planet venus --layout extended
```

Use `--help` to list overrides for the output directory, final time, tolerances, and plot resolution.

## Detailed explanation: how the scripts work

### Execution flow

The files in `scripts/` are deliberately small drivers. They add the local `src/`
directory to Python's import path, select the planet and plot layout, and then call
the common command-line workflow in `pasha_toomre.cli`. The numerical model is not
duplicated between the two scripts.

A complete run proceeds as follows:

1. `config_from_args()` starts with an immutable planet preset and applies any
   command-line overrides. The Venus preset uses `r0 = 0.7226`, `Z0 = -5`, and
   `tau_final = 80`.
2. `initial_state()` builds the eight-component state vector
   `[Z, Vz, x, y, z, vx, vy, vz]`.
3. A short probe integration locates the upward crossing of `Z = 0`. Its physical
   time is stored as `t_c` and defines the encounter-centered time coordinate
   `tau = 6(t - t_c) / pi`.
4. The main integration restarts from the original initial conditions and advances
   through the requested final `tau`. SciPy's adaptive `DOP853` solver supplies a
   dense solution, so reports and plots can use different sampling grids without
   reintegrating the orbit.
5. Scalar event functions detect pericenters and apocenters from zero crossings of
   the three-dimensional radial velocity. Their directions distinguish inward-to-
   outward crossings from outward-to-inward crossings.
6. The dense solution is sampled on an adaptive report grid and on uniform plotting
   grids. The report grid is intentionally much finer near `tau = 0`, where the
   encounter changes most rapidly.
7. The diagnostics module derives orbital quantities, checks numerical invariants,
   and passes named arrays to the reporting and plotting modules.
8. The workflow writes text, CSV, JSON, and PNG products. Any expected validation
   or integration error is reported clearly and causes a nonzero exit status.

The two supplied drivers differ only in presentation:

| Driver | Fixed choices | Result |
|---|---|---|
| `scripts/run_test_venus.py` | Venus, `overview` layout | Original 12 diagnostics in two 6-by-2 figures |
| `scripts/run_venus_extended.py` | Venus, `extended` layout | Adds `eta` and groups 13 diagnostics across four figures |

### Physical model and initial conditions

The Sun and an equal-mass intruder move head-on along the vertical axis. In the
chosen barycentric coordinates, the Sun is at `(0, 0, Z)` and the intruder is at
`(0, 0, -Z)`. The planet is a massless test particle at `(x, y, z)`: it responds to
both stars but does not alter their trajectories.

At the beginning of a Venus run, the planet is on a circular orbit about the Sun:

- its heliocentric position is `(r0, 0, 0)`;
- the Sun begins at `Z0 = -5` and the intruder at `-Z0 = 5`;
- the planet's tangential speed is `sqrt(1 / r0)` in the adopted dimensionless
  units;
- the stellar approach speed is `Vz0 = sqrt(1 / (2 |Z0|))`.

The planet feels ordinary point-mass gravity from both stars. Only the mutual
Sun-intruder force is softened near their exact overlap. The implemented
fourth-order Plummer prescription is

```text
q = Z / R
a_Z = -q (q^2 + 5/8) / [4 R^2 (q^2 + 1/4)^(5/2)]
```

where `R` is the configurable softening scale. This makes the stellar acceleration
finite at maximum approach while recovering the intended large-separation
behavior. The model name is resolved through the softening registry, allowing
additional prescriptions to be added without rewriting the equations of motion or
the solver.

### Encounter time and event handling

The original reconstructed code was blocked by an event callback that returned an
array. SciPy requires a scalar event value. The library now uses a scalar
`sun_crosses_zero` event returning `Z`, marked terminal with positive crossing
direction. This reliably finds the central encounter and prevents ambiguous array
comparisons inside `solve_ivp`.

Pericenter and apocenter detection uses

```text
v_r,3D = (x_rel vx_rel + y_rel vy_rel + z_rel vz_rel) / r
```

with all relative quantities measured from the Sun. A negative-to-positive zero
crossing is a pericenter; a positive-to-negative crossing is an apocenter. Event
times and radii come directly from the adaptive integrator rather than from minima
or maxima on the sampled plot grid.

### Diagnostics

The common diagnostic set contains:

- `Z`: Sun position in the barycentric frame;
- `zeta = z - Z`: the planet's vertical coordinate relative to the Sun;
- `r` and `rxy`: three-dimensional and projected heliocentric radii;
- `phi_phase`: unwrapped orbital phase;
- `vr`, `vphi`, and `vz_rel`: cylindrical radial, tangential, and vertical
  velocities relative to the Sun;
- `vr3d`: the full three-dimensional radial velocity used for extrema events;
- `e` and `ra`: instantaneous osculating eccentricity and aphelion radius;
- `inclination`: inclination derived from the heliocentric angular momentum.

The extended layout inserts `eta = z + Z`, the vertical planet-intruder separation.
Osculating `e` and `ra` are masked before `tau = 0.1` because they are not useful
during the strongest non-Keplerian part of the encounter. Small negative
eccentricity radicands caused only by floating-point roundoff are clamped to zero;
materially negative values raise a diagnostic error instead of being hidden with an
unconditional absolute value.

The code also follows all three angular-momentum components internally. Axial
angular momentum `L_z` should be conserved because the force field remains
axisymmetric. Its initial value, final value, and absolute drift are written to the
JSON summary as an end-to-end accuracy check.

### Sampling and output files

Each run directory contains:

| File | Purpose |
|---|---|
| `trajectory_report.txt` | Human-readable fixed-width table of the sampled diagnostics |
| `trajectory_report.csv` | Machine-readable version suitable for analysis or comparison |
| `extrema.csv` | Integrator-detected pericenters and apocenters with `tau` and radius |
| `run_summary.json` | Configuration, solver statistics, encounter time, and conservation checks |
| `*.png` | Figures for the selected layout |

The text and CSV reports use the same rows. Plot curves are evaluated independently
at higher uniform resolution, so making smoother figures does not enlarge the
tables or change the numerical integration. The overview layout produces two
6-by-2 diagnostic figures. The extended layout produces four grouped figures and
includes the extra `eta` column in both reports.

### Library map and extension points

| Module | Responsibility |
|---|---|
| `config.py` | Validated immutable configuration and Earth/Venus presets |
| `softening.py` | Softened stellar-force implementations and model registry |
| `dynamics.py` | Initial state, equations of motion, and scalar event functions |
| `simulation.py` | Probe and main integrations, time transforms, and sampling grids |
| `diagnostics.py` | Derived orbital quantities, extrema records, and invariant checks |
| `reporting.py` | Text, CSV, extrema, and JSON writers |
| `plotting.py` | Overview and extended Matplotlib layouts |
| `cli.py` | Argument parsing and end-to-end workflow orchestration |
| `errors.py` | Project-specific exceptions converted into clear CLI failures |

To add a planet, define a validated preset in `config.py` or supply the existing
CLI overrides. To add a softening law, implement the same acceleration interface in
`softening.py` and register its name. New diagnostics should be computed once in
`diagnostics.py`; reports and plot layouts can then select them by name. The letter
mentions second-order Plummer and polynomial softening, but their equations are not
present in the supplied documents, so those models must be provided explicitly
rather than inferred.

### Reproducing the supplied Earth plot

The same library can reproduce the Earth calculation used for the JPEG in `docs/`:

```bash
python -m pasha_toomre.cli --planet earth --layout overview
```

This writes to `outputs/earth_overview/`. The reproduced curves and principal
post-encounter extrema agree with the supplied reference plot; the numerical values
are available in its `extrema.csv` and `run_summary.json`.

## Test

```bash
python -m unittest discover -s tests -v
```

All 16 tests pass. They cover the softened force, initial conditions, scalar events, central-encounter integration, diagnostics, output-grid construction, and conservation of axial angular momentum.

## Verified full-run results

Both Venus layouts use identical dynamics and produce identical numerical results through `tau = 80`:

| Quantity | Value |
|---|---:|
| `t_c` at the upward `Z = 0` crossing | 10.540970288572 |
| Initial `tau` | -20.131770316932 |
| Initial `L_z` | 0.850058821494136 |
| Final `L_z` | 0.850058821493976 |
| Absolute `L_z` drift | 1.6065e-13 |
| Main DOP853 function evaluations | 8420 |
| Detected radial extrema | 9 |

The main post-encounter extrema are:

| Type | `tau` | Heliocentric radius |
|---|---:|---:|
| Pericenter | 0.9603687579 | 0.5040131602 |
| Apocenter | 21.5278848504 | 3.9941105088 |
| Pericenter | 41.8366169374 | 0.5135493210 |
| Apocenter | 62.1767821027 | 4.0005199153 |

Complete event lists are stored in each run's `extrema.csv`.

## Outputs

The overview run writes to `outputs/test_venus/`; the extended run writes to `outputs/venus_extended/`. Each directory contains:

- a fixed-width text report;
- a numeric CSV report;
- event-derived extrema in CSV format;
- a JSON run summary;
- the requested PNG figures.

The extended CSV contains 14 columns including `tau`; the overview CSV contains 13. Both contain 329 sampled rows and no unexpected non-finite values. Eccentricity and aphelion are intentionally masked before `tau = 0.1`.

## Softening models

The library implements the fourth-order Plummer law given explicitly in the Venus documents. Its registry allows more models to be added without changing the integrator. The letter mentions second-order Plummer and polynomial prescriptions but does not provide their equations, so the implementation does not invent them.
