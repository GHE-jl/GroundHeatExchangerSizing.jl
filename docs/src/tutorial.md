# Tutorial

## Installation

The GHE-jl packages are not yet registered; develop the ecosystem packages side by side:

```julia
using Pkg
Pkg.develop(path = "../BoreholeResistance.jl")
Pkg.develop(path = "../GroundResponse.jl")
Pkg.develop(path = "../GroundHeatExchanger.jl")
Pkg.develop(path = ".")
Pkg.instantiate()
```

## Sizing from an hourly ground load

Both families expose a dispatcher that takes the **hourly** ground load profile and resamples it
internally for the requested `level` (`:L2`, `:L3` or `:L4`).

```julia
using GroundHeatExchangerSizing

Q  = my_hourly_ground_loads          # 8760 ground loads [W] (– extraction, + rejection)
xy = [0.0 0.0]                       # borehole coordinates (nb × 2) [m]

rb, D          = 0.075, 4.0
ks, Cs         = 2.25, 2.5e6
s, ro, ri      = 0.075, 0.0167, 0.013
kg, Cg         = 1.73, 3.0e6
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

Every sizing call returns a named tuple `(H, H_low, H_high, sol_low, sol_high)`: the governing
length `H` (the larger of the two operating limits), the per-limit lengths, and the two
`Optimization.jl` solution objects.

## Choosing the level and the loads

- `level = :L2` reduces the load to three pulses (yearly / monthly / peak), `:L3` uses 12 monthly
  average and peak loads, `:L4` uses the full 8760-hour profile.
- The peak duration (`tp`, default 6 h) and design period (`ny`, default 10 years) are keywords.
- To call a level directly with pre-resampled loads, use `alternative_sizing_L2` /
  `_L3` / `_L4` (and the `outlet_sizing_*` equivalents); the load-resampling helpers
  ([`Q_hourly_to_three_pulses`](@ref), [`Q_hourly_to_monthly`](@ref), …) produce the required shapes.

## Multiple boreholes

Pass the borehole coordinates as an `nb × 2` matrix; the finite-line-source g-function then accounts
for borehole-to-borehole thermal interaction by spatial superposition. Ground loads `Q` are the
**total** field loads and the flow rate `V` is **per borehole loop**.

```julia
using GroundHeatExchanger: borefield_rectangle
xy = borefield_rectangle(5, 5, 6.0, 6.0)   # 5 × 5 field, 6 m spacing
```
