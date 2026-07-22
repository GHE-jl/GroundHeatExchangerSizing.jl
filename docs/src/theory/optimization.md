# Optimisation

Every sizing equation reduces to the same one-dimensional problem: find the borehole length ``H``
for which the simulated fluid (or outlet) temperature just reaches an operating limit
``T_\text{lim}`` without exceeding it. For each limit the package minimises the residual

```math
\min_{H \in [50, 250]\,\text{m}} \; \bigl| T_\text{lim} - \operatorname{ext} T(H) \bigr|,
```

where the extremum is the **minimum** temperature for a low limit (ground cooling) and the
**maximum** for a high limit (ground heating), split at 10 °C (Dion & Pasquier 2025). For the
three-pulse (L2) equations the temperature is a scalar, so the extremum is the value itself.

The problem is solved with [Optimization.jl](https://github.com/SciML/Optimization.jl) using the
Optim.jl `Fminbox(LBFGS())` backend and `AutoFiniteDiff()` gradients — the same configuration as the
model inversions of ThermalResponseTest.jl. Finite differences are used because the g-functions,
the outlet transfer function and its neural-network short-term half are not dual-number
differentiable.

The search is bounded to ``[50, 250]`` m: the lower bound keeps the design realistic and the upper
bound matches the validity range of the short-term ANN, so both sizing families share identical
bounds. The solver, AD backend, initial guess and bounds are all keyword arguments
(`optimizer`, `adtype`, `H0`, `lb`, `ub`).

A GHE must satisfy both operating limits, so the problem is solved twice and the **governing** design
is the longer of the two lengths, returned as `H` alongside the per-limit lengths `H_low`, `H_high`
and the two `Optimization.jl` solution objects.
