# GroundHeatExchangerSizing.jl

A Julia package to **size vertical ground heat exchangers (GHE)**. It finds the borehole length
required to keep the heat carrier fluid temperature within its operating limits while covering the
ground thermal loads required to maintain setpoint temperature in a building. The package provides two sizing families. Each family works at three levels of
detail (L2 three-pulse, L3 monthly, L4 hourly):

1. **Alternative ASHRAE sizing equation.** This is the g-function form of Ahmadfard & Bernier (2018,
   2019). It removes the temperature-penalty term by evaluating the finite-line-source g-function for
   the actual borefield.
2. **Borehole-outlet transfer-function sizing.** This is the method of Dion & Pasquier (2025). It
   replaces the borehole-wall g-function with a dimensionless transfer function defined at the
   borehole outlet. The transfer function embeds the fluid residence time, the pipe and borehole
   geometry, and the effective borehole resistance.

`GroundHeatExchangerSizing.jl` is a sizing layer of the GHE-jl ecosystem. Its backend is
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl), that re-exports
the [GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) ground models and the
[BoreholeResistance.jl](https://github.com/GHE-jl/BoreholeResistance.jl) resistance and water
property functions. It also provides the temporal-superposition `convolution` and the
`outlet_transfer_function`. The two families solve for the borehole length in different ways. The
alternative ASHRAE equation uses a fixed-point iteration. The borehole-outlet method uses a bounded
one-dimensional optimisation with [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)'s `Brent`
method.

> **Scope.** The package sizes on **ground** thermal loads. If you have only a building load and
> average heating/cooling COPs, `Q_COP(Qb, COP_heating, COP_cooling)` converts them to ground loads
> as a convenience. Load-dependent or temperature-dependent COPs and capacity limiting remain the
> responsibility of `GroundSourceHeatPumpDesign.jl`.

## Quick start

```julia
using GroundHeatExchangerSizing

Q  = my_hourly_ground_loads          # 8760 ground loads [W] (– extraction, + rejection)
xy = [0.0 0.0]                       # borehole coordinates (nb × 2) [m]

# Physical inputs (single U-tube):
rb, D          = 0.075, 4.0          # borehole radius, buried depth [m]
ks, Cs         = 2.25, 2.5e6         # ground conductivity [W/mK], heat capacity [J/m³K]
s, ro, ri      = 0.075, 0.0167, 0.013
kg, Cg         = 1.73, 2.5e6         # grout conductivity / heat capacity
kp, Cp         = 0.40, 1.9e6         # pipe  conductivity / heat capacity (Cp fixed by DeepANN)
kf, cf, ρf, μf = 0.468, 4019.0, 1026.0, 3.37e-3   # fluid properties
V              = 4.9e-4              # flow rate per borehole loop [m³/s]
T0, Tlim       = 10.0, [0.0, 35.0]   # undisturbed ground temp, [low, high] limits [°C]

# Alternative ASHRAE sizing (g-function)
H = alternative_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim; level = :L4)

# Borehole-outlet transfer-function sizing (Dion & Pasquier 2025)
H = outlet_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim; level = :L4)
```

Both families return a single value: the governing borehole length `H` in metres, the larger of the
two operating-limit lengths.

## Sizing equations

Both families expose per-level functions and a dispatcher. The dispatcher takes the **hourly** load
and resamples it internally (`level = :L2 | :L3 | :L4`).

| Function | Level | Load input | Ground response |
|---|---|---|---|
| `alternative_sizing_L2` | three pulses | 3 × 2 pulses | field FLS g-function at 3 superposition times |
| `alternative_sizing_L3` | monthly | 12 × 3 monthly | FLS g-function + convolution |
| `alternative_sizing_L4` | hourly | 8760 | FLS g-function + convolution |
| `alternative_sizing` | dispatcher | hourly | selects the level and resamples |
| `outlet_sizing_L2` | three pulses | 3 × 2 pulses | outlet transfer function at 3 times |
| `outlet_sizing_L3` | monthly | 12 × 3 monthly | outlet transfer function + convolution |
| `outlet_sizing_L4` | hourly | 8760 | outlet transfer function + convolution |
| `outlet_sizing` | dispatcher | hourly | selects the level and resamples |

The finite-line-source g-function is evaluated under the equal-mean-wall-temperature boundary
condition (**BC-II**) by successive spatial superposition, on a `FLSModel`. The effective borehole
thermal resistance `Rb*` is computed with the first-order multipole method and the
axial short-circuit correction (`resistance_ULoop_effective`). It is recomputed at each candidate
length.

### Common arguments

`(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim)`

| Argument | Meaning | Unit |
|---|---|---|
| `Q` | ground load (hourly for the dispatcher; native shape for a level) | W |
| `xy` | borehole coordinates (nb × 2); `[0.0 0.0]` for one borehole | m |
| `rb`, `D` | borehole radius, buried depth | m |
| `ks`, `Cs` | ground thermal conductivity, volumetric heat capacity | W/mK, J/m³K |
| `s`, `ro`, `ri` | shank spacing, pipe outer / inner radius | m |
| `kg`, `Cg` | grout conductivity, volumetric heat capacity (`Cg` for outlet only) | W/mK, J/m³K |
| `kp`, `Cp` | pipe conductivity, volumetric heat capacity (`Cp` for outlet only) | W/mK, J/m³K |
| `kf`, `cf`, `ρf`, `μf` | fluid conductivity, specific heat, density, viscosity | W/mK, J/kgK, kg/m³, kg/m/s |
| `V` | volumetric flow rate in **one** U-tube loop, and **per borehole** | m³/s |
| `T0`, `Tlim` | undisturbed ground temperature, `[low, high]` limits | °C |

The keywords are `level` (`:L2`, `:L3` or `:L4`, dispatcher only), `tp` (peak duration in hours,
default 6), and `ny` (design period in years, default 10).

> **Note.** The borehole-outlet transfer function uses `DeepANN` (Pasquier & Marcotte, 2020), the
> default short-term ANN of `outlet_transfer_function`. That network is valid for `H` between 50
> and 250 m and for a certain interval of geometric and thermal parameters. The length search is
> bounded to the same 50 to 250 m range. Inputs outside the training ranges are clamped by the
> backend, which prints a warning.

## Thermal load analysis

The load-resampling helpers (in `thermal_load_analysis.jl`) convert an hourly profile into the forms
the L2 and L3 equations need. The dispatchers use them internally.

| Function | Purpose |
|---|---|
| `Q_analysis(Q)` | Summarise an hourly profile into a `QLoads` structure (yearly / monthly / peak) |
| `Q_hourly_to_monthly(Q)` | Hourly → 12 × 3 monthly average / cooling-peak / heating-peak loads |
| `Q_hourly_to_three_pulses(Q)` | Hourly → 3 × 2 three-pulse loads (yearly, monthly, peak; cooling & heating) |
| `Q_monthly_to_hourly(Qm, tp)` | Monthly loads → 8760 × 2 hourly profile (peak over the last `tp` seconds) |
| `Q_monthly_to_three_pulses(Qm)` | Monthly loads → 3 × 2 three-pulse loads |
| `Q_cutoff(Q, ch, cc)` | Scale the heating/cooling peaks of an hourly profile |

## Solving for the length

The two families reach the borehole length in different ways.

The alternative ASHRAE equation uses a fixed-point iteration. It starts from an initial length and
computes a new length from the closed-form ASHRAE expression. It repeats this until the length stops
changing by more than 0.01 m or a maximum iteration count is reached. The g-function is recomputed
at each iterate because it depends on the length.

The borehole-outlet method uses a bounded optimisation. For each operating limit it minimises the
absolute difference between the limit and the temperature extremum over `H` between 50 and 250 m.
The residual is a smooth, one-dimensional, unimodal function of `H`, so the solver is Optim.jl's
`Brent` method. Brent is derivative-free. It brackets the minimum inside the bounds and shrinks the
bracket until the length stops moving by more than 1 mm. It needs no initial guess and no gradient,
so it converges in a few tens of evaluations.

Both families size against a low and a high operating limit. The governing length is the larger of
the two per-limit results.

## Scripts

Runnable, plotted examples live in `script/` (its own environment). Run them from the package root:

```
julia --project=script -e 'using Pkg; Pkg.instantiate()'
julia --project=script script/script_alternative_sizing.jl
```

Each script displays its figure.

| Script | What it shows |
|---|---|
| `script_alternative_sizing.jl` | Alternative ASHRAE L2/L3/L4 lengths on the four Ahmadfard & Bernier cases |
| `script_outlet_sizing.jl` | Outlet transfer-function sizing vs the alternative equation, per level |
| `script_thermal_load_analysis.jl` | The hourly, monthly and three-pulse load conversions |
| `Ahmadfard_cases.jl` | The four reference cases (loads and parameters) of Ahmadfard & Bernier (2019) |

## Installation

The GHE-jl packages are not yet registered. Clone the ecosystem packages side by side and develop
them locally:

```julia
using Pkg
Pkg.develop(path = "../BoreholeResistance.jl")
Pkg.develop(path = "../GroundResponse.jl")
Pkg.develop(path = "../GroundHeatExchanger.jl")
Pkg.develop(path = ".")
Pkg.instantiate()
```

## Dependencies

### Library

| Package | Used for |
|---|---|
| [GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) | FLS g-functions and outlet transfer function, temporal-superposition `convolution` (re-exports GroundResponse.jl and BoreholeResistance.jl) |
| [GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) | `FLSModel`, `ground_response`, spatial superposition |
| [BoreholeResistance.jl](https://github.com/GHE-jl/BoreholeResistance.jl) | Effective borehole resistance `Rb*`, water properties |
| [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) | Bounded one-dimensional `Brent` optimisation of the borehole length (outlet method) |

### Scripts only

| Package | Used in |
|---|---|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | All visualisation scripts |
| [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) | `script_alternative_sizing.jl`, `script_outlet_sizing.jl` timing blocks |

## References

- Ahmadfard, M., & Bernier, M. (2018). Modifications to ASHRAE's sizing method for vertical ground
  heat exchangers. *Science and Technology for the Built Environment*, 24(7), 803–817.
  https://doi.org/10.1080/23744731.2018.1423816
- Ahmadfard, M., & Bernier, M. (2019). A review of vertical ground heat exchanger sizing tools
  including an inter-model comparison. *Renewable and Sustainable Energy Reviews*, 110, 247–265.
  https://doi.org/10.1016/j.rser.2019.04.045
- Dion, G., & Pasquier, P. (2025). Ground heat exchanger sizing using borehole outlet transfer
  function. *Science and Technology for the Built Environment*.
  https://doi.org/10.1080/23744731.2025.2523200
- Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of artificial neural networks to
  near-instant construction of short-term g-functions. *Applied Thermal Engineering*, 143, 910–921.
  https://doi.org/10.1016/j.applthermaleng.2018.07.137
- Pasquier, P., & Marcotte, D. (2020). Robust identification of volumetric heat capacity and
  analysis of thermal response tests by Bayesian inference with correlated residuals. *Applied
  Energy*, 261, 114394. https://doi.org/10.1016/j.apenergy.2019.114394
