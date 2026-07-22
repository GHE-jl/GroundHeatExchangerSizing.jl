# Overview

Sizing a vertical ground heat exchanger means finding the borehole length ``H`` for which the
heat-pump fluid temperature reaches an operating limit ``T_\text{lim}`` over the design period
without exceeding it. Following Spitler & Bernier (2016) and Ahmadfard & Bernier (2019), sizing tools
are classified by the temporal detail of the ground load:

1. **L2** uses three pulses (the yearly average, the monthly average during the peak month, and the
   peak). They are applied over 10 years, one month and 4 to 6 hours.
2. **L3** uses monthly average and monthly peak loads (36 values per year).
3. **L4** uses the full 8760-hour load profile.

All levels rest on the same superposition of ground loads with a ground thermal response. This
package provides two choices of that response, both backed by
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl):

1. the **finite-line-source g-function** at the borehole wall. See the
   [Alternative ASHRAE equation](@ref).
2. the dimensionless **borehole-outlet transfer function**. See the
   [Outlet transfer function](@ref) method.

In both cases the borehole length satisfies a temperature-limit condition. The alternative equation
finds it by a fixed-point iteration. The outlet method finds it by a bounded
[Optimisation](@ref).

## Sign convention and inputs

Ground loads ``Q`` are negative for heat extraction (ground cooling, building heating) and positive
for heat rejection (ground heating, building cooling). A GHE must satisfy a lower operating limit
(dominated by the cooling loads) and a higher limit (dominated by the heating loads). The governing
design is the longer of the two lengths.

The package sizes on **ground** loads. Converting building loads to ground loads through the
heat-pump coefficient of performance is the responsibility of `GroundSourceHeatPumpDesign.jl`.
