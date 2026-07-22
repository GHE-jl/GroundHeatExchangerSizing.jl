module GroundHeatExchangerSizing

using GroundHeatExchanger

# Source files (included once, in dependency order)
include("thermal_load_analysis.jl")     # hourly ↔ monthly ↔ three-pulse ground-load resampling
include("utils.jl")                     # deviation metrics
include("alternative_sizing.jl")        # alternative ASHRAE sizing equation (fixed-point iteration)
include("outlet_sizing.jl")             # borehole-outlet transfer-function sizing (optimisation)

# Thermal-load analysis and resampling
export QLoads,
    Q_analysis,
    Q_hourly_to_monthly,
    Q_hourly_to_three_pulses,
    Q_monthly_to_hourly,
    Q_monthly_to_three_pulses,
    Q_cutoff

# Alternative ASHRAE sizing equation (Ahmadfard & Bernier 2018, 2019)
export alternative_sizing,
    alternative_sizing_L2,
    alternative_sizing_L3,
    alternative_sizing_L4

# Borehole-outlet transfer-function sizing (Dion & Pasquier 2025)
export outlet_sizing,
    outlet_sizing_L2,
    outlet_sizing_L3,
    outlet_sizing_L4

# Utilities
export Q_COP, deviation_metrics

end
