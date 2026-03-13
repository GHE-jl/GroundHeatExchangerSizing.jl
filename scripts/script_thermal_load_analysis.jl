"""
Script that showcases the functions of `QAnalysis.jl` to analyze the thermal load profiles of a GHE sizing.
The script is structured in four main sections:
1. `Q_hourly_to_monthly`: Resample the hourly heat load profile to monthly average and monthly peak loads
2. `Q_hourly_to_three_pulse`: Resample the hourly heat load profile to three pulse loads
3. `Q_monthly_to_hourly`: Resample the monthly average and monthly peak loads to hourly values for a year.
4. `Q_monthly_to_three_pulse`: Resample the monthly average and monthly peak loads to three pulse loads.
"""

using CairoMakie

includet("../src/GHESizings.jl")
using .GHESizings

# 1. Select case
includet("Ahmadfard_cases.jl")
Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = ahmadfard_cases(4)

# 2. Analyse the heat load profile
Q = Q_analysis(Q₀)

# 3. Test other functions for conversion
Q_h_m = Q_hourly_to_monthly(Q₀)
Q_h_3p = Q_hourly_to_three_pulse(Q₀)
Q_m_h = Q_monthly_to_hourly(Q_h_m, 6 * 3600)
Q_m_3p = Q_monthly_to_three_pulse(Q_h_m)

