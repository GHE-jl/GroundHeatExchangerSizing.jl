```@meta
CurrentModule = GroundHeatExchangerSizing
```

# GroundHeatExchangerSizing.jl

Sizing of vertical **ground heat exchangers (GHE)**: finding the borehole length that keeps the
heat-pump fluid temperature within its operating limits while covering the ground thermal loads.

Two sizing families are provided, each at three levels of detail (L2 three-pulse, L3 monthly,
L4 hourly):

- **Alternative ASHRAE sizing equation** — the g-function form of Ahmadfard & Bernier (2018, 2019),
  which removes the temperature-penalty term by evaluating the finite-line-source g-function for the
  actual borefield.
- **Borehole-outlet transfer-function sizing** — Dion & Pasquier (2025), which replaces the
  borehole-wall g-function by a dimensionless transfer function defined at the borehole outlet.

The package is a sizing layer of the GHE-jl ecosystem. Its backend is
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) (which re-exports
[GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) and
[BoreholeResistance.jl](https://github.com/GHE-jl/BoreholeResistance.jl)), and the borehole length is
found with [Optimization.jl](https://github.com/SciML/Optimization.jl).

The package sizes on **ground** thermal loads; building-to-ground load conversion (heat-pump COP)
belongs to `GroundSourceHeatPumpDesign.jl`.

See the [Tutorial](@ref) to get started, the theory section for the governing equations, and the
[API reference](@ref) for every function.
