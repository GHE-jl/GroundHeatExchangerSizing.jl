# Outlet transfer function

Dion & Pasquier (2025) replace the borehole-**wall** g-function with a dimensionless transfer
function ``\bar g`` defined at the borehole **outlet**. The outlet is the entering water temperature
on the heat-pump source side. The outlet temperature is the temporal superposition of the ground
loads with this transfer function (Eq. 7 of Dion & Pasquier 2025):

```math
T_\text{out}(t) = T_g + \sum_i \frac{Q_i - Q_{i-1}}{\dot V\,C_f}\,\bar g(t - t_{i-1}, H)
```

Here ``\dot V`` is the volumetric flow rate and ``C_f`` the fluid volumetric heat capacity. Unlike
the g-function form, there is **no explicit borehole resistance term**. The effect of ``R_b^\ast`` is
carried inside ``\bar g``.

## Building the transfer function

``\bar g`` is the `outlet_transfer_function` of GroundHeatExchanger.jl. The short-term half is the
artificial neural network of Pasquier et al. (2018). The long-term half is the FLS borehole-wall
response converted to a dimensionless outlet transfer function through the effective borehole
resistance:

```math
\bar g = \left(\frac{g_\text{wall}}{2\pi k_s} + R_b^\ast\right)\frac{\dot V\,C_f}{H}
```

The two halves are spliced at the 7-day contact time. The ANN embeds the fluid residence time, the
pipe and borehole geometry, and the borehole thermal capacity. Its geometric input is the **shank
spacing** ``s``, which is distinct from the buried depth ``D`` used by the long-term `FLSModel`.

## The three levels

**L2 (three pulses).** Evaluate ``\bar g`` at the three superposition times and form
``T_\text{out} = T_g + (Q_y\Gamma_y + Q_m\Gamma_m + Q_h\Gamma_h)/(N_b\dot V C_f)`` with
``\Gamma_y = \bar g(t_3) - \bar g(t_2)``, ``\Gamma_m = \bar g(t_2) - \bar g(t_1)`` and
``\Gamma_h = \bar g(t_1)``.

**L3 and L4.** Superimpose the monthly-expanded or hourly load with ``\bar g`` by FFT `convolution`.
The full L4 profile is solved in the spectral domain without load aggregation.

The transfer function is per borehole (field-average) and ``\dot V`` is per borehole, so the
borehole-count factor cancels. Total ground loads ``Q`` are used with ``Q/(N_b\dot V C_f)``.

## Validity range

The short-term ANN is trained for ``H \in [110, 200]`` m and for a narrow band of geometric and
thermal parameters (see `outlet_transfer_function`). The length search is bounded to the same
``[110, 200]`` m range. Inputs outside the training ranges are clamped by the backend, which prints a
warning.
