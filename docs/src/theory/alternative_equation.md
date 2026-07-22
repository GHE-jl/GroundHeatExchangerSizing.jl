# Alternative ASHRAE equation

The alternative ASHRAE sizing equation (Ahmadfard & Bernier 2018, 2019) writes the mean fluid
temperature as the temporal superposition of the ground loads with a finite-line-source (FLS)
g-function evaluated for the actual borefield. The field g-function already contains the
borehole-to-borehole thermal interaction, so the temperature-penalty term of the classical ASHRAE
equation is not needed:

```math
T_f(t) = T_g + \frac{1}{N_b H}\left[\sum_i Q_i\,R_{g,i} + Q_h\,R_b^\ast\right]
```

Here ``T_g`` is the undisturbed ground temperature, ``N_b`` the number of boreholes, ``H`` the
borehole length, ``Q_i`` the ground load pulses [W], ``R_{g,i}`` the effective ground thermal
resistances built from the field g-function [°C·m/W], and ``R_b^\ast`` the effective borehole
thermal resistance.

## Ground thermal resistances

The FLS g-function ``g(t, H)`` [°C·m/W] is evaluated under the equal-mean-wall-temperature boundary
condition (**BC-II**) by successive spatial superposition, on a `FLSModel(H, D, k_s, C_s)`.

**L2 (three pulses).** The superposition times are ``t_1 = t_h``, ``t_2 = t_m + t_h`` and
``t_3 = t_y + t_m + t_h`` (peak ``t_h``, month ``t_m``, design period ``t_y``). The resistances are:

```math
R_{g,h} = g(t_1), \qquad
R_{g,m} = g(t_2) - g(t_1), \qquad
R_{g,y} = g(t_3) - g(t_2)
```

**L3 and L4.** The load is superimposed with the full g-function by FFT `convolution`. For L3 the
monthly load is first expanded to hourly, with the peak over the final ``t_p`` hours of each month.

```math
T_f(t) = T_g + \big(q \ast g\big)(t) + q(t)\,R_b^\ast,
\qquad q(t) = \frac{Q(t)}{N_b H}
```

## Effective borehole resistance

``R_b^\ast`` is computed with the first-order multipole method and the axial thermal short-circuit
correction (`resistance_ULoop_effective`, single U-tube). It depends on ``H`` through the fluid
thermal-capacity term, so it is recomputed at every candidate length.

## Solving for the length

For each operating limit the length satisfies ``T_f = T_\text{lim}`` (the peak fluid temperature for
L2, the extremum of the hourly series for L3 and L4). The equation is solved by a fixed-point
iteration. The iteration starts from an initial length and computes a new length from the closed-form
expression above. It repeats until the length stops changing by more than 0.01 m or a maximum
iteration count is reached. The g-function is recomputed at each iterate because it depends on the
length. The vector of iterates is returned as `Hi`. The governing design is the longer of the low
and high limit lengths.
