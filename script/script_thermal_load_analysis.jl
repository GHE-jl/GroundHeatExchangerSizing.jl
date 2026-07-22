# Thermal load analysis and resampling. Shows the conversions used to feed the sizing equations:
#   1. Q_hourly_to_monthly       — hourly → 12 monthly average / cooling-peak / heating-peak loads
#   2. Q_hourly_to_three_pulses  — hourly → yearly / monthly / peak three-pulse loads
#   3. Q_monthly_to_hourly       — monthly loads → hourly profile (peak over the last tp hours)
#   4. Q_monthly_to_three_pulses — monthly loads → three-pulse loads
#
# Run from the package root with the script environment:
#   julia --project=script script/script_thermal_load_analysis.jl

using GroundHeatExchangerSizing
using CairoMakie

include("Ahmadfard_cases.jl")

p = ahmadfard_cases(4)
Q = Float64.(p.Q)

analysis = Q_analysis(Q)                          # QLoads summary structure
Qm = Q_hourly_to_monthly(Q)                       # 12 × 3 [Qma Qmc Qmh]
Q3 = Q_hourly_to_three_pulses(Q)                  # 3 × 2 [Qy Qy; Qm_c Qm_h; Qh_c Qh_h]
Qhr = Q_monthly_to_hourly(Qm, 6 * 3600.0)         # 8760 × 2 [cooling, heating]
Q3m = Q_monthly_to_three_pulses(Qm)               # 3 × 2

println("Yearly average load : $(round(analysis.y; digits = 1)) W")
println("Peak heating / cooling : $(round(analysis.hₕ; digits = 1)) / ",
    "$(round(analysis.hₗ; digits = 1)) W")

fig = Figure(size = (820, 520))
ax1 = Axis(fig[1, 1], title = "Hourly ground load", xlabel = "Hour", ylabel = "Q [W]")
lines!(ax1, 1:length(Q), Q)
ax2 = Axis(fig[1, 2], title = "Monthly loads", xlabel = "Month", ylabel = "Q [W]")
barplot!(ax2, 1:12, Qm[:, 1]; label = "average")
scatter!(ax2, 1:12, Qm[:, 2]; label = "cooling peak")
scatter!(ax2, 1:12, Qm[:, 3]; label = "heating peak")
axislegend(ax2; position = :rb)
ax3 = Axis(fig[2, 1:2], title = "Monthly-to-hourly reconstruction (heating column)",
    xlabel = "Hour", ylabel = "Q [W]")
lines!(ax3, 1:size(Qhr, 1), Qhr[:, 2])

mkpath(joinpath(@__DIR__, "figures"))
save(joinpath(@__DIR__, "figures", "thermal_load_analysis.png"), fig)
println("Saved figures/thermal_load_analysis.png")
