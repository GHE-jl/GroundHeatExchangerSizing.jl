# Optimisation

The borehole-outlet method finds the length by a bounded optimisation. The goal is to find the
borehole length ``H`` for which the simulated outlet temperature just reaches an operating limit
``T_\text{lim}`` without exceeding it. For each limit the package minimises the residual

```math
\min_{H \in [110, 200]\,\text{m}} \; \bigl| T_\text{lim} - \operatorname{ext} T(H) \bigr|.
```

The extremum is the **minimum** temperature for a low limit (ground cooling) and the **maximum** for
a high limit (ground heating), split at 10 °C (Dion & Pasquier 2025). For the three-pulse (L2)
equations the temperature is a scalar, so the extremum is the value itself.

The problem is solved with [Optimization.jl](https://github.com/SciML/Optimization.jl) using the
Optim.jl `Fminbox(LBFGS())` backend and `AutoFiniteDiff()` gradients. Finite differences are used
because the g-functions, the outlet transfer function and its neural-network short-term half are not
dual-number differentiable.

The search is bounded to ``[110, 200]`` m. This is the training range of the short-term ANN in the
borehole length. Outside this range the ANN clamps ``H``, so the short-term transfer function and the
residual become flat in ``H`` and the optimiser cannot make progress. Keeping the search inside the
ANN range keeps the residual smooth and avoids that stall. The bounds are fixed by design and are not
keyword arguments.

A GHE must satisfy both operating limits, so the problem is solved twice. The **governing** design is
the longer of the two lengths, returned as `H` alongside the per-limit lengths `H_low` and `H_high`.

The alternative ASHRAE equation does not use this optimiser. It solves for the length by a
fixed-point iteration instead. See the [Alternative ASHRAE equation](@ref).
