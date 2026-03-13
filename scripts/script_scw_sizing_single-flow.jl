"""
Script performing SCW sizing using an adaptation of the alternative sizing equation for constant
flow rates.
"""

# using BenchmarkTools
using CairoMakie

includet("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")
includet("../src/GHESizings.jl")
using .GroundHeatExchanger, .GHESizings

# Select case
includet("Ahmadfard_cases.jl")
Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = ahmadfard_cases(4)

# Ensure the yearly heat load profile is in the right type and scale it
Q = Float64.(Q₀) * 0.19

# Define other variable to simulate the SCW
k = (s = 2.5,           # Ground thermal conductivity
    p = 0.4)            # Pipe thermal conductivity
Cs = 2.4e6              # Ground volumetric heat capacity
K = 1e-7                # Assumption: Groundwater hydraulic conductivity
r = (b = 0.076,         # Borehole radius
    o = 0.03,           # Pipe outlet radius
    i = 0.025)          # Pipe inlet radius
V = 100.0/6e4           # Circulating flow rate in the SCW (100 L/min)
β = 0.0                 # Bleed ratio (between 0 and 1)
Cf = water_cp(T.g) * water_ρ(T.g) # Fluid volumetric heat capacity at T.g

# Sizing parameters
tp = 6.0                # Time of the peak heat pulse (typically either 4 or 6 hours)
ny = 10                 # Number of years to use in sizing
Tlim = [0.0, 35.0]      # Ground operating temperature limit

# Three pulses SCW sizing equation with constant flow rates
Q_3p = Q_hourly_to_three_pulses(Q)
@time Ĥ, out = scw_sizing_three_pulses(Q_3p, Tlim, tp, k.s, k.p, Cs, K, r.b, r.o, r.i, V, β, T.g)

# Monthly SCW sizing equation with constant flow rates
Q_m = Q_hourly_to_monthly(Q)
@time Ĥ, out = scw_sizing_monthly(Q_m, Tlim, ny, tp, k.s, k.p, Cs, K, r.b, r.o, r.i, V, β, T.g)

# Hourly SCW sizing equation with constant flow rates
@time Ĥ, out = scw_sizing_hourly(Q, Tlim, ny, k.s, k.p, Cs, K, r.b, r.o, r.i, V, β, T.g)

# Compute the optimized SCW response
t = 3600.0:3600.0:(ny * 365 * 24 * 3600.0)
g = βils_outlet(t, [Ĥ k.s Cs], k.p, [Ĥ K], r.b, r.o, r.i, Ĥ, Ĥ-10.0, V, β, T.g)  # Output in W/m

# Compute the inlet and outlet temperature response
Q_ = repeat(Q, ny) # Repeat the yearly profile for n years
T_out = T.g .+ convolution(impulse_func(Q_) / Ĥ, g)
T_in = T_out .+ Q_ / (V * Cf)

# Figures
fig = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))

ax = Axis(fig[1, 1], xlabel = L"$t$ (d)", ylabel = L"$Q$ (W)")
lines!(ax, t[1:365*24] / (3600 * 24), Q, linewidth = 1.5)

ax = Axis(fig[1, 2], xlabel = L"$t$ (y)", ylabel = L"$g$ (°Cm/W)", xscale = log10)
lines!(ax, t / (3600 * 24 * 365), g, linewidth = 1.5)

ax = Axis(fig[2, 1:2], xlabel = L"$t$ (y)", ylabel = L"$T$ (°C)")
lines!(ax, t / (3600 * 24 * 365), T_in, linewidth = 1.5, label = "Inlet Temp.")
lines!(ax, t / (3600 * 24 * 365), T_out, linewidth = 1.5, label = "Outlet Temp.")
hlines!(ax, Tlim, color = :black, linestyle = :dash, linewidth = 1.5)
axislegend(ax; position = :rb)

display(fig);