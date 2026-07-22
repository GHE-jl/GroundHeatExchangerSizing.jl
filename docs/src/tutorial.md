# Tutorial

## Installation

The GHE-jl packages are not yet registered. Develop the ecosystem packages side by side:

```julia
using Pkg
Pkg.develop(path = "../BoreholeResistance.jl")
Pkg.develop(path = "../GroundResponse.jl")
Pkg.develop(path = "../GroundHeatExchanger.jl")
Pkg.develop(path = ".")
Pkg.instantiate()
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
kp, Cp         = 0.40, 1.54e6
kf, cf, ρf, μf = 0.468, 4019.0, 1026.0, 3.37e-3
V              = 4.9e-4
T0, Tlim       = 10.0, [0.0, 35.0]

# Alternative ASHRAE sizing (g-function)
res = alternative_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim;
                         level = :L4)

# Borehole-outlet transfer-function sizing (needs the grout/pipe heat capacities Cg, Cp)
res = outlet_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim;
                    level = :L4)

res.H          # governing borehole length [m]
```

The two families return different named tuples. `alternative_sizing` returns `(H, Hi)`, where `H` is
the governing length (the larger of the two operating limits) and `Hi` is the vector of fixed-point
iterates. `outlet_sizing` returns `(H, H_low, H_high)`, where `H` is the governing length and
`H_low`, `H_high` are the two per-limit lengths.

## Choosing the level and the loads

The `level` keyword sets how the load is reduced. `:L2` reduces the load to three pulses (yearly,
monthly and peak). `:L3` uses 12 monthly average and peak loads. `:L4` uses the full 8760-hour
profile.

The peak duration (`tp`, default 6 h) and the design period (`ny`, default 10 years) are keywords.

You can also call a level directly with pre-resampled loads. Use `alternative_sizing_L2`, `_L3` or
`_L4` (and the `outlet_sizing_*` equivalents). The load-resampling helpers
([`Q_hourly_to_three_pulses`](@ref), [`Q_hourly_to_monthly`](@ref) and others) produce the required
shapes.

## Multiple boreholes

Pass the borehole coordinates as an `nb × 2` matrix. The finite-line-source g-function then accounts
for borehole-to-borehole thermal interaction by spatial superposition. The ground loads `Q` are the
**total** field loads and the flow rate `V` is **per borehole loop**.

```julia
using GroundHeatExchanger: borefield_rectangle
xy = borefield_rectangle(5, 5, 6.0, 6.0)   # 5 × 5 field, 6 m spacing
```
