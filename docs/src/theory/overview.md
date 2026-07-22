# Overview

Sizing a vertical ground heat exchanger means finding the borehole length ``H`` for which the
heat-pump fluid temperature reaches — but does not exceed — an operating limit ``T_\text{lim}`` over
the design period. Following Spitler & Bernier (2016) and Ahmadfard & Bernier (2019), sizing tools
are classified by the temporal detail of the ground load:

- **L2** — three pulses (yearly average, monthly average during the peak month, and the peak),
  applied over 10 years, one month and 4–6 hours;
- **L3** — monthly average and monthly peak loads (36 values per year);
- **L4** — the full 8760-hour load profile.

All levels rest on the same superposition of ground loads with a ground thermal response. This
package provides two choices of that response, both backed by
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl):

1. the **finite-line-source g-function** at the borehole wall — the
   [Alternative ASHRAE equation](@ref);
2. the dimensionless **borehole-outlet transfer function** — the
   [Outlet transfer function](@ref) method.

In both cases the borehole length is the root of a temperature-limit condition, found by the shared
bounded [Optimisation](@ref).

## Sign convention and inputs

Ground loads ``Q`` are negative for heat extraction (ground cooling, building heating) and positive
for heat rejection (ground heating, building cooling). A GHE must satisfy a lower operating limit
(dominated by the cooling loads) and a higher limit (dominated by the heating loads); the governing
design is the longer of the two lengths.

The package sizes on **ground** loads. Converting building loads to ground loads through the
heat-pump coefficient of performance is the responsibility of `GroundSourceHeatPumpDesign.jl`.
