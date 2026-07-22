```@meta
CurrentModule = GroundHeatExchangerSizing
```

# GroundHeatExchangerSizing.jl

Sizing of vertical **ground heat exchangers (GHE)**. The goal is to find the borehole length that
keeps the heat-pump fluid temperature within its operating limits while covering the ground thermal
loads.

The package provides two sizing families. Each family works at three levels of detail (L2
three-pulse, L3 monthly, L4 hourly):

1. **Alternative ASHRAE sizing equation.** This is the g-function form of Ahmadfard & Bernier (2018,
   2019). It removes the temperature-penalty term by evaluating the finite-line-source g-function for
   the actual borefield.
2. **Borehole-outlet transfer-function sizing.** This is the method of Dion & Pasquier (2025). It
   replaces the borehole-wall g-function with a dimensionless transfer function defined at the
   borehole outlet.

The package is a sizing layer of the GHE-jl ecosystem. Its backend is
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl), which re-exports
[GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) and
[BoreholeResistance.jl](https://github.com/GHE-jl/BoreholeResistance.jl). The alternative ASHRAE
equation solves for the length by a fixed-point iteration. The borehole-outlet method solves for the
length with [Optimization.jl](https://github.com/SciML/Optimization.jl).

The package sizes on **ground** thermal loads. Converting building loads to ground loads through the
heat-pump COP belongs to `GroundSourceHeatPumpDesign.jl`.

See the [Tutorial](@ref) to get started, the theory section for the governing equations, and the
[API reference](@ref) for every function.
