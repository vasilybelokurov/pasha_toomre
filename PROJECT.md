# Sun–Intruder Encounter and Planetary-Orbit Response

## Project overview

This project studies how a close, axial encounter between the Sun and a second massive body (the “intruder”) perturbs a planetary orbit. The immediate numerical problem is to regularize the Sun–intruder force during their maximum approach, where an unsoftened point-mass interaction would become singular. Several gravitational-softening prescriptions are to be tested and compared so that a numerically stable model can be chosen without significantly distorting the resulting planetary orbit.

The current code is a Venus adaptation of an earlier Earth calculation. It follows the Sun and a massless test planet in a symmetric Sun–intruder geometry, integrates the encounter with a high-order adaptive solver, records a detailed orbital time series, locates radial extrema, checks conservation of the axial angular momentum, and produces overview and close-encounter plots.

The supplied material consists of two Word documents containing versions of the Python program, a letter to Alan Toomre describing the Earth tests and their physical interpretation, and a JPEG example of the diagnostic plots.

## Scientific objective

The main objective is to select an appropriate model for smoothing the Sun–intruder interaction near their nominal collision at the origin. A successful prescription should:

1. remove the force singularity at maximum approach;
2. allow stable and accurate numerical integration through the central encounter;
3. preserve the correct large-separation gravitational force;
4. minimize artificial changes in the planet’s post-encounter orbital elements;
5. give converged results under tighter integration tolerances and finer sampling;
6. preserve quantities that should remain conserved by axial symmetry, especially the planet’s $L_z$.

The longer-term goal, stated in the letter to Alan Toomre, is to extend the calculation from individual Earth and Venus tests to all planets and produce comparative plots of their orbital responses.

## Dynamical model

### Geometry

The calculation uses a non-rotating Cartesian frame centered on the Sun–intruder center of mass. The motion is axially symmetric about the $z$-axis:

- the Sun is at $(0,0,Z)$;
- the intruder is at $(0,0,-Z)$;
- the planet is at $(x,y,z)$.

The Sun begins below the origin, so $Z_0<0$, and moves upward. The intruder follows the mirror trajectory. At $Z=0$, the two massive bodies reach their nominal point of maximum approach. The code calls this the “collapse” event and defines it as $\tau=0$.

The planet is treated as a massless test particle. It feels the gravity of both massive bodies but does not affect either of them. The Sun and intruder have equal normalized gravitational masses in the planet’s equations of motion.

### State vector

The integrated state is

\[
\boldsymbol{y}=(Z,V_Z,x,y,z,v_x,v_y,v_z).
\]

The first two components describe the Sun. The remaining six describe the planet in the center-of-mass frame.

### Units

The code uses the standard normalization $G M_\odot=1$, with distances expressed in astronomical units. This interpretation is supported by the adopted values:

- $r_0=0.7226$, the initial Venus orbital radius in AU;
- $R=0.00464913034$, approximately the solar radius in AU.

In these units, the time unit is $\sqrt{\mathrm{AU}^3/(G M_\odot)}=1/(2\pi)$ years. The plotted time coordinate is therefore

\[
\tau=\frac{6}{\pi}(t-t_c),
\]

where $t_c$ is the moment at which $Z=0$. Consequently, one unit of $\tau$ is one Earth month, matching the convention cited in the letter and in the earlier LT76 work.

## Venus initial conditions

The supplied programs use

\[
Z_0=-5, \qquad r_0=0.7226, \qquad R=0.00464913034.
\]

The initial upward speed of the Sun is

\[
V_{Z,0}=\sqrt{\frac{1}{2|Z_0|}},
\]

which is the parabolic, zero-energy value for the adopted symmetric infall model at large separation.

The planet initially lies in the Sun’s horizontal plane and shares its vertical motion:

\[
(x_0,y_0,z_0)=(r_0,0,Z_0),
\]

\[
(v_{x,0},v_{y,0},v_{z,0})=\left(0,\sqrt{\frac{1}{r_0}},V_{Z,0}\right).
\]

Thus, relative to the Sun, Venus initially has the circular speed appropriate to radius $r_0$.

## Gravitational forces

### Softened Sun–intruder force

The current Venus code labels its regularization as a fourth-order Plummer force. Defining

\[
q=\frac{Z}{R},
\]

the Sun’s acceleration is

\[
\ddot Z=-\frac{q\left(q^2+5/8\right)}
{4R^2\left(q^2+1/4\right)^{5/2}}.
\]

This force is finite and odd in $Z$, so it passes smoothly through zero rather than diverging. At large separation, $|Z|\gg R$, it approaches

\[
\ddot Z\simeq-\frac{\operatorname{sgn}(Z)}{4Z^2},
\]

which is the expected acceleration of either massive body when their separation is $2|Z|$. Near the origin it becomes linear, $\ddot Z\simeq-5Z/R^3$, and therefore removes the point-mass singularity.

The Alan Toomre letter describes earlier Earth experiments with two other regularizations: a second-order Plummer model and a smooth polynomial force inside $|Z|<R$. Those Earth runs gave very similar results; the quoted difference in final aphelion radius was about $0.0016$.

### Planetary acceleration

Define the planet’s displacement from the Sun and intruder as

\[
\boldsymbol{r}_{\rm S}=(x,y,z-Z), \qquad
\boldsymbol{r}_{\rm I}=(x,y,z+Z),
\]

with magnitudes $r_{\rm S}$ and $r_{\rm I}$. The planet obeys

\[
\ddot{\boldsymbol{r}}=-\frac{\boldsymbol{r}_{\rm S}}{r_{\rm S}^3}
-\frac{\boldsymbol{r}_{\rm I}}{r_{\rm I}^3}.
\]

Only the mutual Sun–intruder motion is softened in the supplied model. The planet’s direct forces from the Sun and intruder remain point-mass forces and would still be singular in an exact planet–body collision.

## Numerical method

The intended program performs two integrations with SciPy’s `solve_ivp` and its explicit eighth-order DOP853 method.

### Probe integration

The first pass searches for the physical crossing time $t_c$ at which the Sun moves upward through $Z=0$. Its nominal tolerances are

- relative tolerance: $10^{-11}$;
- absolute tolerance: $10^{-13}$.

The detected $t_c$ establishes the shifted monthly time coordinate $\tau$.

### Main integration

The second pass runs from the supplied initial state through

\[
t_{\rm end}=t_c+\frac{\pi}{6}\tau_{\rm fin},
\qquad \tau_{\rm fin}=80,
\]

using tighter tolerances:

- relative tolerance: $10^{-12}$;
- absolute tolerance: $10^{-14}$.

Dense output is enabled so that the state can be evaluated at arbitrary requested values of $\tau$.

### Radial-extremum events

The heliocentric three-dimensional radial velocity is

\[
v_{r,3D}=\frac{xv_x+yv_y+(z-Z)(v_z-V_Z)}{r_{\rm S}}.
\]

Zeros crossed from negative to positive identify pericenters, while zeros crossed from positive to negative identify apocenters. The event solver is intended to determine their times and radii more precisely than a search on the output grid.

### Output sampling

The report grid is deliberately adaptive in $\tau$. It is coarse far from the encounter and progressively refined near $\tau=0$, reaching a nominal spacing of $0.0004$ in the central interval. This dense sampling is separate from DOP853’s own adaptive internal steps.

## Reported quantities

The intended text report and plots contain the following diagnostics:

- $Z$: Sun’s vertical coordinate in the center-of-mass frame;
- $\zeta=z-Z$: planet’s vertical displacement from the Sun;
- $\eta=z+Z$: planet’s vertical displacement from the intruder, present only in one code variant;
- $r=r_{\rm S}$: heliocentric three-dimensional distance;
- $r_{xy}=\sqrt{x^2+y^2}$: horizontal distance from the symmetry axis;
- $\phi_{\rm phase}=[\operatorname{atan2}(y,x)/(2\pi)]\bmod 1$: azimuthal phase in turns;
- $v_r=(xv_x+yv_y)/r_{xy}$: radial velocity in the horizontal plane;
- $v_\phi=(xv_y-yv_x)/r_{xy}$: azimuthal velocity;
- $v_{z,\rm rel}=v_z-V_Z$: planet’s vertical velocity relative to the Sun;
- $v_{r,3D}$: heliocentric orbital radial velocity;
- $e$: eccentricity of the instantaneous heliocentric osculating orbit;
- $r_a$: aphelion radius of that osculating orbit when $e<1$;
- $i$: inclination of the osculating orbit in degrees.

The angular momentum is calculated from the heliocentric relative position and velocity. In particular,

\[
L_z=xv_y-yv_x,
\]

and

\[
i=\operatorname{atan2}\!\left(\sqrt{L_x^2+L_y^2},L_z\right).
\]

The osculating Kepler quantities use

\[
e=\sqrt{1+\left(v_{\rm rel}^2-\frac{2}{r}\right)L^2},
\qquad
r_a=\frac{L^2}{1-e}\quad(e<1).
\]

The current code applies an absolute value inside the square root and suppresses $e$ and $r_a$ before $\tau=0.1$. These values must be interpreted as instantaneous heliocentric osculating elements. They are not conserved while the intruder’s perturbation remains strong, especially during the first few positive months.

## Conservation check

Because both massive bodies remain on the $z$-axis, their forces exert no torque about that axis. The planet’s $L_z$ should therefore remain constant even though the total angular-momentum magnitude and the orbital inclination can change substantially.

The program compares the initial and final $L_z$. The letter reports agreement through the first 12 decimal places in the Earth runs. This is a valuable internal check, but it does not by itself establish the physical fidelity of a softening model: different central force laws can all preserve axial symmetry while producing slightly different energies and post-encounter orbital elements.

## Physical interpretation from the Earth tests

The letter to Alan Toomre gives the following preliminary interpretation of the Earth calculations:

1. By $\tau\approx-1$, the inclination has already reached roughly $5^\circ$. It then rises rapidly toward $90^\circ$ as the Sun moves upward relative to the planet near the central encounter.
2. A small kink appears in the heliocentric distance $r$ near $\tau=0$, but not in $r_{xy}$. The same feature appeared with both polynomial and Plummer softening, so it is unlikely to arise solely from derivative discontinuities in the polynomial prescription.
3. After the crossing, the planet remains close to the symmetry axis for more than a month. During this interval it gains substantial orbital energy and horizontal angular momentum.
4. The letter interprets this phase as being driven primarily by the intruder: the planet–intruder vertical separation remains small, and the planet temporarily resembles an object moving on a nearly circular horizontal orbit around the intruder.
5. At about $\tau\approx2$, the planet separates from the intruder and returns to a predominantly heliocentric orbit. Its eccentricity, inclination, and aphelion radius then settle toward their post-encounter values.
6. The final Earth orbit shows inclination oscillations of approximately $27^\circ$–$29^\circ$.

For the two Earth softening prescriptions, the letter quotes the following heliocentric radial extrema. Parenthesized values refer to polynomial softening:

| Extremum | Plummer result | Polynomial result |
|---|---:|---:|
| First pericenter | $r=0.697076$, $\tau=1.56244$ | $r=0.697594$, $\tau=1.56330$ |
| First apocenter | $r=5.513360$, $\tau=34.92605$ | $r=5.511777$, $\tau=34.91888$ |
| Second pericenter | $r=0.710225$, $\tau=67.87211$ | $r=0.710450$, $\tau=67.85401$ |

These numbers describe the Earth tests, not a verified Venus run.

## Supplied files and code variants

### `docs/to AT.docx`

This is the English letter to Alan Toomre. It documents the Earth experiments, the comparison between second-order Plummer and polynomial softening, the numerical extrema, and the preliminary physical interpretation summarized above.

### `docs/IMG_3091BA26C461-1.jpeg`

This image is a 12-panel, full-range diagnostic plot over $\tau\in[-10,80]$, with the central crossing marked at $\tau=0$. Its format matches the two-canvas, 6-by-2 plotting design. In context, it appears to be reference output associated with the Earth work rather than a verified output from the supplied Venus documents.

### `docs/2026 07 27 Test Venus PLUM. tau80, Z0 = - 5.0.docx`

This variant contains the 12 standard plotted quantities and intends to generate two 6-by-2 figures:

1. the full interval $[-10,80]$;
2. a detailed central interval $[-2,4]$.

It saves the intended figures as `trajectory_full_overview.png` and `trajectory_singularity_core.png`.

### `docs/2026 07 27 Venus PLUM. tau80, Z0 = - 5.0.docx`

This variant introduces the additional quantity

\[
\eta=z+Z,
\]

the planet’s vertical displacement from the intruder. It intends to divide 13 diagnostics across four figures: two for the full interval and two for the central interval.

The feature description in the colleague’s message—an extra column and reformatted plots—matches this variant more closely than the file named `Test Venus`. Attachment order cannot be established from the filenames alone. Moreover, $\eta$ was added to the plotting arrays and labels but not to the text-report table, so the requested extra table column is not actually complete.

## Present condition of the Python code

Neither Word document currently contains a directly executable Python program. The physical and numerical intent is recoverable, and most of the integration/reporting section is shared between the files, but several defects must be repaired during extraction to `.py` files.

### Confirmed blocking defects

1. The zero-crossing event returns the entire state vector:

   ```python
   def sun_crosses_zero(t, state):
       return state
   ```

   A SciPy event function must return a scalar. For the intended event, this must return `state[0]`, the Sun’s $Z$-coordinate. As written, the probe integration cannot locate $t_c$ correctly and should fail when SciPy evaluates the event sign.

2. The code after the radial-extremum loop has lost essential Python structure in both documents. Indentation is missing, multiple statements are joined with Unicode line-separator characters, and prose headings or separator lines lack comment markers.

3. Several exponent expressions in the plotting function have been converted into forms such as `x2`, `rxy2`, `Lx2`, and `vx2` instead of `x**2`, `rxy**2`, `Lx**2`, and `vx**2`.

4. One subtraction uses a typographic en dash (`z – Z`) instead of the Python minus sign (`z - Z`).

5. In the `Test Venus` plotting section, labels such as `Figure 1:` and raw separator lines occur as executable text, causing syntax errors. The title also contains the malformed literal `TAU_finВ`, including a Cyrillic character, rather than the numerical upper limit or a correctly formatted variable.

6. In the 13-diagnostic variant, the four-figure plotting block has collapsed into one paragraph beginning with `#`. If copied literally, most or all of that block becomes a single comment instead of executable code.

7. The four-figure variant calls methods such as `axes2.legend(...)` on an array of axes and calls `fig2.delaxes(axes2)` or `fig4.delaxes(axes4)` with an axes array. The intended operations require a specific axes object; for example, the unused eighth panel would be `axes2[7]` or `axes4[7]`.

8. The additional $\eta$ diagnostic is absent from the written report’s header and rows, despite the request for an extra tabular column.

### Numerical and interpretive cautions

- The probe crossing event should normally be terminal once the upward $Z=0$ crossing is found. Leaving it non-terminal and integrating hundreds of time units beyond the required event wastes work and could complicate event selection if later crossings occur.
- Applying `abs` inside the eccentricity square root can conceal an invalid negative radicand caused by numerical or formula errors. Small round-off excursions should be handled explicitly and large negative values should be treated as failures.
- The chosen softening modifies only the Sun–intruder orbit. Convergence tests should separately establish that the planet never approaches either point mass closely enough for its unsoftened force to dominate the numerical error.
- A conserved $L_z$ tests axial symmetry and integration accuracy, but comparison of softening prescriptions should also monitor energy transfer, pericenter/apocenter times and radii, post-encounter $e$, $r_a$, $i$, and convergence with tolerance and softening scale.

## Recommended reconstruction workflow

The safest continuation is:

1. extract the common physical model into a real Python module rather than executing code copied from Word;
2. repair the scalar crossing event and reconstruct all lost operators, line breaks, comments, and indentation;
3. implement the $\eta$ quantity consistently in both the report table and plotting data;
4. retain the requested four-figure formatting if that is the intended “second” variant;
5. add explicit checks for integration success, event detection, finite values, and valid output intervals;
6. reproduce the known Earth reference case before accepting the Venus results;
7. verify $L_z$ conservation and compare the recovered Earth extrema with the values in the letter;
8. run the Venus case and produce both full-range and central-range figures;
9. compare alternative softening laws on a controlled parameter grid using identical tolerances and initial conditions;
10. choose the preferred softening model based on convergence and stability, not only visual similarity of the plots.

## Project status

The two Word-embedded programs have been reconstructed as a shared Python library with two runnable Venus drivers. The implementation repairs the event, indentation, operator, table, and subplot defects described above. It also adds the missing intruder-relative `eta` column to the extended report.

Both full runs now complete through `tau=80`. Sixteen automated tests pass, the two layouts produce identical dynamics, and the measured axial angular-momentum drift is approximately `1.61e-13`. Reports, summaries, extrema, and six visually inspected plots are available under `outputs/`. See `README.md` for verified commands and numerical results.
