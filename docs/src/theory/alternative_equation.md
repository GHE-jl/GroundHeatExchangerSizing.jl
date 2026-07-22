# Alternative ASHRAE equation

The alternative ASHRAE sizing equation (Ahmadfard & Bernier 2018, 2019) writes the mean fluid
temperature as the temporal superposition of the ground loads with a finite-line-source (FLS)
g-function evaluated for the actual borefield. Because the field g-function already contains the
borehole-to-borehole thermal interaction, the temperature-penalty term of the classical ASHRAE
equation is not needed:

```math
T_f(t) = T_g + \frac{1}{N_b H}\left[\sum_i Q_i\,R_{g,i} + Q_h\,R_b^\ast\right]
```

where ``T_g`` is the undisturbed ground temperature, ``N_b`` the number of boreholes, ``H`` the
borehole length, ``Q_i`` the ground load pulses [W], ``R_{g,i}`` the effective ground thermal
resistances built from the field g-function [°C·m/W], and ``R_b^\ast`` the effective borehole
thermal resistance.

## Ground thermal resistances

The FLS g-function ``g(t, H)`` [°C·m/W] is evaluated under the equal-mean-wall-temperature boundary
condition (**BC-II**) by successive spatial superposition, on a `FLSModel(H, D, k_s, C_s)`.

**L2 (three pulses).** With superposition times ``t_1 = t_h``, ``t_2 = t_m + t_h`` and
``t_3 = t_y + t_m + t_h`` (peak ``t_h``, month ``t_m``, design period ``t_y``):

```math
R_{g,h} = g(t_1), \qquad
R_{g,m} = g(t_2) - g(t_1), \qquad
R_{g,y} = g(t_3) - g(t_2)
```

**L3 / L4.** The monthly (expanded to hourly, with the peak over the final ``t_p`` hours of each
month) or hourly load is superimposed with the full g-function by FFT `convolution`:

```math
T_f(t) = T_g + \big(q \ast g\big)(t) + q(t)\,R_b^\ast,
\qquad q(t) = \frac{Q(t)}{N_b H}
```

## Effective borehole resistance

``R_b^\ast`` is computed with the first-order multipole method and the axial thermal short-circuit
correction (`resistance_ULoop_effective`, single U-tube). It depends on ``H`` through the fluid
thermal-capacity term and is therefore recomputed at every candidate length.

## Solving for the length

For each operating limit the length is the root of ``T_f = T_\text{lim}`` (the peak fluid
temperature for L2, the extremum of the hourly series for L3/L4), obtained by the bounded
[Optimisation](@ref). The governing design is the longer of the low- and high-limit lengths.
