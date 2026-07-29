# Sun–Intruder Simulation Library Design

Date: 2026-07-28
Status: Implemented and verified

## Purpose

Reconstruct the two Python programs embedded in the supplied Word documents as a reusable, tested simulation library. The library will reproduce both documented Venus output layouts, remove the known execution blockers, and provide a stable base for future planets and gravitational-softening models.

## Goals

- Preserve the documented Sun–intruder–planet dynamical model and Venus initial conditions.
- Put the common physics, integration, diagnostics, reporting, and plotting code in one package.
- Provide two runnable drivers corresponding to the two Word-document variants.
- Add the missing intruder-relative diagnostic, `eta = z + Z`, consistently to the extended report and plots.
- Run both full Venus simulations through `tau = 80`.
- Verify numerical success, event detection, finite outputs, extrema, and axial angular-momentum conservation.
- Produce deterministic reports and publication-quality PNG plots without requiring an interactive display.
- Make new planets and softening laws configurable without duplicating the solver.

## Non-goals

- Do not invent the polynomial softening law mentioned in the Alan Toomre letter; its formula is not supplied.
- Do not claim that the reconstructed Venus results reproduce the quoted Earth extrema.
- Do not include planetary back-reaction on the Sun or intruder.
- Do not soften the planet’s direct force from either massive body unless a later model explicitly requests it.
- Do not build a general-purpose N-body framework.

## Package structure

```text
pyproject.toml
src/pasha_toomre/
    __init__.py
    config.py
    softening.py
    dynamics.py
    simulation.py
    diagnostics.py
    reporting.py
    plotting.py
    cli.py
scripts/
    run_test_venus.py
    run_venus_extended.py
tests/
    test_softening.py
    test_dynamics.py
    test_simulation.py
    test_diagnostics.py
```

The package will depend only on NumPy, SciPy, and Matplotlib. Tests will use the standard-library `unittest` framework so the project needs no additional test dependency.

## Configuration

A frozen `SimulationConfig` dataclass will contain:

- run name and output directory;
- planet name and initial heliocentric radius;
- initial Sun coordinate `Z0`;
- solar smoothing scale `R`;
- final monthly time `tau_final`;
- solver method and probe/main tolerances;
- output-grid and plotting bounds;
- softening-law identifier;
- diagnostics layout;
- whether to display figures interactively.

Named constructors will provide the documented Venus configuration and an Earth configuration for later reference runs. Drivers may override individual values through command-line arguments.

## Softening interface

`softening.py` will expose a registry mapping a name to a scalar acceleration function with signature:

```python
acceleration(Z: float | ndarray, R: float) -> float | ndarray
```

The initial implementation will include `plummer4`, reconstructed from the documents:

```text
q = Z / R
a_Z = -q (q^2 + 5/8) / [4 R^2 (q^2 + 1/4)^(5/2)]
```

Tests will verify that this acceleration:

- is finite at the origin;
- is odd in `Z`;
- points toward the origin;
- approaches `-sign(Z)/(4 Z^2)` at large separation.

The registry will provide a clear extension point for the second-order Plummer and polynomial prescriptions once their exact formulas are available.

## Dynamics

The integrated state remains:

```text
[Z, Vz, x, y, z, vx, vy, vz]
```

The Sun uses the selected softened acceleration. The planet uses the sum of the unsoftened point-mass accelerations from the Sun at `(0, 0, Z)` and the intruder at `(0, 0, -Z)`.

The dynamics layer will validate state shape and reject non-finite values or zero planet–body separations with explicit errors.

Event functions will be scalar:

- `sun_crosses_zero` returns `state[0]` and detects an upward crossing;
- `find_pericenter` returns heliocentric `v_r_3D` with direction `+1`;
- `find_apocenter` returns heliocentric `v_r_3D` with direction `-1`.

This removes the shared blocking error in both Word documents.

## Integration flow

1. Build the documented initial state. Venus begins at `(r0, 0, Z0)` with velocity `(0, sqrt(1/r0), Vz0)`, where `Vz0 = sqrt[1/(2|Z0|)]`.
2. Run a probe DOP853 integration with the zero-crossing event marked terminal.
3. Check `sol.success`, require exactly one usable upward crossing, and record `t_c`.
4. Compute `t_end = t_c + tau_final pi/6`.
5. Run the main DOP853 integration with dense output and pericenter/apocenter events.
6. Check solver success, solution bounds, finite state values, and event arrays.
7. Generate the documented adaptive report grid, with its finest sampling around `tau = 0`.
8. Evaluate diagnostics, extrema, and conservation metrics.
9. Write outputs and create plots.

The plotting grids `[-10, tau_final]` and `[-2, 4]` will be validated against the integrated solution interval before dense-output evaluation.

## Diagnostics

The library will calculate:

- `Z`;
- heliocentric vertical separation `zeta = z - Z`;
- intruder-relative vertical separation `eta = z + Z`;
- heliocentric distance `r`;
- horizontal distance `r_xy`;
- azimuthal phase;
- horizontal radial and azimuthal velocities;
- vertical velocity relative to the Sun;
- heliocentric three-dimensional radial velocity;
- angular-momentum components and magnitude;
- inclination;
- heliocentric osculating eccentricity and aphelion.

The base layout omits `eta` from its report and plots to match `Test Venus`. The extended layout includes it in both.

For eccentricity, the code will replace the original unconditional `abs` with an explicit tolerance:

- clamp a radicand in `[-epsilon, 0)` to zero;
- raise a diagnostic error for a substantially negative radicand.

As in the source documents, eccentricity and aphelion will be masked before `tau = 0.1`. The report will identify them as heliocentric osculating quantities.

## Outputs

Each driver will write to its own directory under `outputs/`.

Common files:

- `trajectory_report.txt`: human-readable fixed-width report;
- `trajectory_report.csv`: machine-readable diagnostic table;
- `extrema.csv`: event-derived pericenters and apocenters;
- `run_summary.json`: configuration, crossing time, solver statistics, and conservation metrics.

The base driver will produce:

- `trajectory_full_overview.png`;
- `trajectory_singularity_core.png`.

The extended driver will produce:

- `trajectory_full_coordinates.png`;
- `trajectory_full_velocities.png`;
- `trajectory_detailed_coordinates.png`;
- `trajectory_detailed_velocities.png`.

Matplotlib will use the non-interactive `Agg` backend by default. A `--show` option will permit interactive display.

## Command-line interface

The general CLI will support commands equivalent to:

```bash
python -m pasha_toomre.cli --planet venus --layout overview
python -m pasha_toomre.cli --planet venus --layout extended
```

The two scripts will be thin, directly runnable wrappers around these configurations:

```bash
python scripts/run_test_venus.py
python scripts/run_venus_extended.py
```

Useful options will include `--output-dir`, `--tau-final`, `--rtol`, `--atol`, and `--show`. Drivers will retain the documented defaults.

## Error handling

The library will raise clear, typed errors for:

- an unknown softening model;
- invalid configuration values;
- failure to detect the Sun’s upward zero crossing;
- unsuccessful SciPy integration;
- evaluation outside the integrated interval;
- non-finite states or diagnostics;
- a planet–body distance below the numerical collision threshold;
- a materially invalid eccentricity radicand.

The CLI will catch these errors, print a concise message, and exit with a nonzero status. Successful runs will print the output directory and key summary values.

## Validation and debugging

Automated tests will cover:

1. the fourth-order Plummer force’s symmetry, regularity, direction, and asymptote;
2. the initial Venus state and circular heliocentric velocity;
3. scalar event functions and event directions;
4. a shortened smoke integration through the central encounter;
5. successful crossing detection near `Z = 0`;
6. diagnostic identities, including `eta = z + Z` and `L_z = x vy - y vx`;
7. finite report data and consistent column lengths;
8. correct plot-data shapes for 12- and 13-diagnostic layouts;
9. preservation of `L_z` within a conservative numerical threshold.

After the tests pass, both full Venus drivers will run through `tau = 80`. The generated summaries and plots will be inspected. The final report will record:

- crossing time;
- solver success and evaluation counts;
- initial/final `L_z` and their difference;
- all detected extrema;
- output file paths;
- any remaining scientific cautions.

## Acceptance criteria

The work is complete when:

- both drivers exit successfully in the available project Python environment;
- all automated tests pass;
- both runs detect the upward `Z = 0` crossing;
- all expected reports and PNG files exist and are non-empty;
- the base report has the documented 12 diagnostics;
- the extended report and plots include `eta`;
- no NaN or infinity occurs outside intentionally masked osculating quantities;
- the reported `L_z` drift is within the tested tolerance;
- every generated plot is visually inspected for missing data, malformed labels, overlap, or unused panels;
- the original Word files remain unchanged.
