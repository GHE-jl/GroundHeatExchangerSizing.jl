# GroundHeatExchangerSizing.jl

A Julia package to **size vertical ground heat exchangers (GHE)** — to find the borehole length
required to keep the heat-pump fluid temperature within its operating limits while covering the
ground thermal loads. Two complementary sizing families are provided, each at three levels of
detail (L2 three-pulse, L3 monthly, L4 hourly):

1. **Alternative ASHRAE sizing equation** — the g-function form of Ahmadfard & Bernier (2018, 2019),
   which removes the temperature-penalty term by evaluating the finite-line-source g-function for
   the actual borefield.
2. **Borehole-outlet transfer-function sizing** — Dion & Pasquier (2025), which replaces the
   borehole-wall g-function with a dimensionless transfer function defined at the borehole outlet,
   embedding the fluid residence time, pipe/borehole geometry and effective borehole resistance.

`GroundHeatExchangerSizing.jl` is a sizing layer of the GHE-jl ecosystem: its backend is
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl), which re-exports the
[GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) ground models and the
[BoreholeResistance.jl](https://github.com/GHE-jl/BoreholeResistance.jl) resistance and water-
property functions, and provides the temporal-superposition `convolution` and the
`outlet_transfer_function`. The borehole length is found with
[Optimization.jl](https://github.com/SciML/Optimization.jl) (Optim.jl `Fminbox(LBFGS())` backend).

> **Scope** — the package sizes on **ground** thermal loads. If you have only a building load and
> average heating/cooling COPs, `Q_COP(Qb, COP_heating, COP_cooling)` converts them to ground loads
> as a convenience. Load- or temperature-dependent COPs and capacity limiting remain the
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
kg, Cg         = 1.73, 3.0e6         # grout conductivity / heat capacity
kp, Cp         = 0.40, 1.54e6        # pipe  conductivity / heat capacity
kf, cf, ρf, μf = 0.468, 4019.0, 1026.0, 3.37e-3   # fluid properties
V              = 4.9e-4              # flow rate per borehole loop [m³/s]
T0, Tlim       = 10.0, [0.0, 35.0]   # undisturbed ground temp, [low, high] limits [°C]

# --- Alternative ASHRAE sizing (g-function) ---
res = alternative_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim;
                         level = :L4)
res.H          # governing borehole length [m]

# --- Borehole-outlet transfer-function sizing (Dion & Pasquier 2025) ---
res = outlet_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim;
                    level = :L4)
```

Every sizing call returns `(H, H_low, H_high, sol_low, sol_high)`: the governing length `H` [m]
(the larger of the two limits), the per-limit lengths, and the two `Optimization.jl` solutions.

## Sizing equations

Both families expose per-level functions and a dispatcher that takes the **hourly** load and
resamples it internally (`level = :L2 | :L3 | :L4`).

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
thermal resistance `Rb*` is always computed internally with the first-order multipole method and the
axial short-circuit correction (`resistance_ULoop_effective`), recomputed at each candidate length.

### Common arguments

`(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, [Cg,] kp, [Cp,] kf, cf, ρf, μf, V, T0, Tlim; kwargs...)`

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
| `V` | volumetric flow rate in one U-tube loop (per borehole) | m³/s |
| `T0`, `Tlim` | undisturbed ground temperature, `[low, high]` limits | °C |

Keywords: `level` (`:L2`/`:L3`/`:L4`, dispatcher only), `tp` (peak duration [h], default 6),
`ny` (design period [years], default 10), `model` (`Rb*` boundary condition `"UHF"`/`"UBW"`/`"mean"`),
and any of `H0`, `lb`, `ub`, `optimizer`, `adtype` forwarded to the optimiser.

> **Note** The borehole-outlet transfer function uses the short-term ANN of Pasquier et al. (2018),
> valid for `H ∈ [110, 200] m` and a narrow band of geometric/thermal parameters; the length search
> is bounded to `[50, 250] m`. Out-of-range inputs are clamped by the backend (with a warning).

## Thermal load analysis

The load-resampling helpers (in `thermal_load_analysis.jl`) convert an hourly profile into the forms
the L2/L3 equations need, and are used internally by the dispatchers.

| Function | Purpose |
|---|---|
| `Q_analysis(Q)` | Summarise an hourly profile into a `QLoads` structure (yearly / monthly / peak) |
| `Q_hourly_to_monthly(Q)` | Hourly → 12 × 3 monthly average / cooling-peak / heating-peak loads |
| `Q_hourly_to_three_pulses(Q)` | Hourly → 3 × 2 three-pulse loads (yearly, monthly, peak; cooling & heating) |
| `Q_monthly_to_hourly(Qm, tp)` | Monthly loads → 8760 × 2 hourly profile (peak over the last `tp` seconds) |
| `Q_monthly_to_three_pulses(Qm)` | Monthly loads → 3 × 2 three-pulse loads |
| `Q_cutoff(Q, ch, cc)` | Scale the heating/cooling peaks of an hourly profile |

## Optimisation

Each sizing method builds, for both temperature limits, a residual `|Tlim − extremum(T(H))|` and
minimises it over `H ∈ [50, 250] m`. The default solver is the Optim.jl `Fminbox(LBFGS())` with
`AutoFiniteDiff()` gradients (the g-functions, transfer functions and neural network are not
dual-number differentiable). The governing length is the larger of the two per-limit results.

## Scripts

Runnable, plotted examples live in `script/` (its own environment). Run them from the package root:

```
julia --project=script -e 'using Pkg; Pkg.instantiate()'
julia --project=script script/script_alternative_sizing.jl
```

Each script saves its figure to `script/figures/`.

| Script | What it shows |
|---|---|
| `script_alternative_sizing.jl` | Alternative ASHRAE L2/L3/L4 lengths on the four Ahmadfard & Bernier cases |
| `script_outlet_sizing.jl` | Outlet transfer-function sizing vs the alternative equation, per level |
| `script_thermal_load_analysis.jl` | The hourly ↔ monthly ↔ three-pulse load conversions |
| `Ahmadfard_cases.jl` | The four reference cases (loads + parameters) of Ahmadfard & Bernier (2019) |

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
| [Optimization.jl](https://github.com/SciML/Optimization.jl) / [OptimizationOptimJL.jl](https://github.com/SciML/Optimization.jl) | Bounded borehole-length optimisation |
| [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl) | Finite-difference gradients for the optimiser |

### Scripts only

| Package | Used in |
|---|---|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | All visualisation scripts |

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
