module GroundHeatExchangerSizing

using Revise

# Functions to convert thermal loads between hourly, monthly, and three-pulse for use in sizing.
includet("thermal_load_analysis.jl")
# Functions to perform the alternative ASHRAE sizing method.
includet("alternative_sizing.jl")
# Functions to perform the borehole outlet sizing method.
includet("borehole_outlet_sizing.jl")
# Functions to perform the SCW sizing method.
includet("standing_column_well_sizing.jl")

# From thermal_load_analysis.jl
export QLoads,
    Q_analysis,
    Q_hourly_to_monthly,
    Q_hourly_to_three_pulses,
    Q_monthly_to_hourly,
    Q_monthly_to_three_pulses,
    Q_cutoff,
    Q_COP,
    heat_pump_performance

# From heat_pump.jl
export heat_pump_performance

# From AlternativeSizing
export alternative_sizing,
    alternative_sizing_L2,
    alternative_sizing_L3,
    alternative_sizing_L4

# From OutletSizing.jl
export outlet_sizing_hourly

# From SCWSizing.jl
export scw_sizing_three_pulses,
    scw_sizing_monthly,
    scw_sizing_hourly
end