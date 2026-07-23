# Optimisation

The borehole-outlet method finds the length by a bounded optimisation. The goal is to find the
borehole length ``H`` for which the simulated outlet temperature just reaches an operating limit
``T_\text{lim}`` without exceeding it. For each limit the package minimises the residual

```math
\min_{H \in [50, 250]\,\text{m}} \; \bigl| T_\text{lim} - \operatorname{extremum} T(H) \bigr|.
```

The extremum is the **minimum** temperature for a low limit (ground cooling) and the **maximum** for
a high limit (ground heating). For the three-pulse (L2)
equations the temperature is a scalar, so the extremum is the value itself.

The residual is a smooth, one-dimensional, unimodal function of ``H``, so the problem is solved with
Optim.jl's `Brent` method. Brent is a derivative-free bounded solver. It brackets the minimum inside
the search window and shrinks the bracket by parabolic interpolation, with a golden-section
fallback, until the length stops moving by more than the tolerance (here 1 mm). It needs no initial
guess and no gradient, which suits this problem because the g-functions, the outlet transfer
function and its neural-network short-term half are not dual-number differentiable. A one-dimensional
bracketing method also converges in a few tens of evaluations, so it avoids the very large
evaluation counts a bound-constrained gradient method spends when the minimum lies on a bound.

The search is bounded to ``[50, 250]`` m. This is the training range of `DeepANN` in the borehole
length (Pasquier & Marcotte 2020) — the default short-term ANN `outlet_transfer_function` evaluates.
Outside this range the ANN clamps ``H``, so the short-term transfer function and the residual become
flat in ``H`` and the optimiser cannot make progress. Keeping the search inside the ANN range keeps
the residual smooth and avoids that stall. The bounds are fixed by design and are not keyword
arguments.

A GHE must satisfy both operating limits, so the problem is solved twice. The **governing** design
is the longer of the two lengths, and only that value `H` is returned.

The alternative ASHRAE equation does not use this optimiser. It solves for the length by a
fixed-point iteration instead. See the [Alternative ASHRAE equation](@ref).
