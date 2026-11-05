# Tutorial

## Installation

`GroundHeatExchangerSizing.jl` is registered in the General registry, along with its dependencies
`GroundHeatExchanger.jl`, `GroundResponse.jl` and `BoreholeResistance.jl`, so they all resolve
automatically:

```julia
using Pkg
Pkg.add("GroundHeatExchangerSizing")
```

## Sizing from an hourly ground load

Both families expose a dispatcher. The dispatcher takes the **hourly** ground load profile and
resamples it internally for the requested `level` (`:L2`, `:L3` or `:L4`).

```julia
using GroundHeatExchangerSizing

Q  = my_hourly_ground_loads          # 8760 ground loads [W] (– extraction, + rejection)
xy = [0.0 0.0]                       # borehole coordinates (nb × 2) [m]

rb, D          = 0.075, 4.0
ks, Cs         = 2.25, 2.5e6
s, ro, ri      = 0.075, 0.0167, 0.013
kg, Cg         = 1.73, 2.5e6
kp, Cp         = 0.40, 1.9e6         # Cp fixed by DeepANN
kf, cf, ρf, μf = 0.468, 4019.0, 1026.0, 3.37e-3
V              = 4.9e-4
T0, Tlim       = 10.0, [0.0, 35.0]

# Alternative ASHRAE sizing (g-function); Tlim at the ground inlet (heat-pump LFT) by default
H = alternative_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim; level = :L4)

# Borehole-outlet transfer-function sizing (needs the grout/pipe heat capacities Cg, Cp)
H = outlet_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim; level = :L4)
```

Both families return a single value: the governing borehole length `H` in metres (the larger of the
two operating-limit lengths).

## Choosing the level and the loads

The `level` keyword sets how the load is reduced. `:L2` reduces the load to three pulses (yearly,
monthly and peak). `:L3` uses 12 monthly average and peak loads. `:L4` uses the full 8760-hour
profile.

The peak duration (`tp`, default 6 h) and the design period (`ny`, default 10 years) are keywords.

For both `alternative_sizing` and `outlet_sizing`, `Tlim_at` sets where `Tlim` applies, seen from
the ground (not the heat pump): `:inlet` (ground inlet = heat-pump LFT, default, the most
conservative), `:outlet` (ground outlet = heat-pump EFT) or `:mean` (mean fluid temperature).

You can also call a level directly with pre-resampled loads. Use `alternative_sizing_L2`, `_L3` or
`_L4` (and the `outlet_sizing_*` equivalents). The load-resampling helpers
([`Q_hourly_to_three_pulses`](@ref), [`Q_hourly_to_monthly`](@ref) and others) produce the required
shapes.

## Multiple boreholes

Pass the borehole coordinates as an `nb × 2` matrix. The finite-line-source g-function then accounts
for borehole-to-borehole thermal interaction by spatial superposition. The ground loads `Q` are the
**total** field loads and the flow rate `V` is **per borehole loop**.

```julia
using GroundHeatExchanger: borefield
xy = borefield(:rectangle, 5, 5, 6.0, 6.0)   # 5 × 5 field, 6 m spacing
```
