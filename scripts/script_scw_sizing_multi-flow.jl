"""
Script performing SCW sizing using an adaptation of the alternative sizing equation for multiple
flow rates.
"""

# using BenchmarkTools
using CairoMakie

includet("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")
includet("../src/GHESizings.jl")
using .GroundHeatExchanger, .GHESizings

# Select case
includet("Ahmadfard_cases.jl")
Q₀, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, T, ~, ~ = ahmadfard_cases(4)

# Ensure the yearly heat load profile is in the right type and scale it
# Q₀ = Float64.(Q₀) * 0.25
# Q = Q_analysis(Q₀)      # Outputs the heat load analysis for three-pulse and monthly sizing
Q₁ = 15000 * ones(365 * 24) # For testing purposes, use a constant heat load profile

Q = Q₁

# Define other variable to simulate the SCW
k = (s = 2.5,           # Ground thermal conductivity
    p = 0.4)            # Pipe thermal conductivity
Cs = 2.4e6              # Ground volumetric heat capacity
K = 1e-7                # Assumption: Groundwater hydraulic conductivity
r = (b = 0.076,         # Borehole radius
    o = 0.03,           # Pipe outlet radius
    i = 0.025)          # Pipe inlet radius
V = [75.0, 125.0, 175.0] / 6e4 # Circulating flow rates in the SCW
β = [0.0, 0.05, 0.10]    # Bleed ratio (between 0 and 1)

# Sizing parameters
tp = 6.0                # Time of the peak heat pulse (typically either 4 or 6 hours)
ny = 10                 # Number of years to use in sizing
Tlim = [0.0, 35.0]      # Ground operating temperature limit
t = 3600.0:3600.0:(ny * 365 * 24 * 3600.0)

# Hourly SCW sizing equation with constant flow rates
T_out = []
T_in = []
Ĥ, out = scw_sizing_hourly(Q, Tlim, ny, k.s, k.p, Cs, K, r.b, r.o, r.i, V, β, T.g)

# Compute the SCW transfer function using the βILS model in GHEModels.jl with the optimal length
t_ = t[set_nodes(length(t), 100)] # Subsample time vector for interpolation
g = zeros(length(t), length(V) + length(β) - 1)  # Preallocate
for i in eachindex(V)
    if i != length(V)
        tmp = βils_outlet(t_, [Ĥ k.s Cs], k.p, [Ĥ K], r.b, r.o, r.i, Ĥ, Ĥ-10.0, V[i], β[1], T.g)
        g[:, i] = pchip_interpolation(t_, tmp, t) # From GHEModels.jl (Utils)
    else
        for j in eachindex(β)
            tmp = βils_outlet(t_, [Ĥ k.s Cs], k.p, [Ĥ K], r.b, r.o, r.i, Ĥ, Ĥ-10.0, V[i], β[j], T.g)
            g[:, i + j - 1] = pchip_interpolation(t_, tmp, t) # From GHEModels.jl (Utils)
        end
    end
end

# Compute temperature response of the SCW using the states optimization algorithm
state = ones(Int, length(t))                # Initial state (flow rate index)
T_out = similar(t)
T_in  = similar(t)
maxH = argmax([out[1].minimizer, out[2].minimizer]) # Dominant phase obtained after sizing

aux_energy, state, T_in, T_out, Q_aux = states_optim!(T_in, T_out, state, repeat(Q, ny) / Ĥ, g, V,
    T.g; mode = Tlim[maxH] < 17.5 ? :heating : :cooling, ΔT = 0.3, maxiter = 25)

# Obtain circulating and bleed flow rates from the state vector
V_circ, V_bleed = zeros(length(t)), zeros(length(t))
for i in eachindex(state)
    if state[i] <= length(V)
        V_circ[i] = V[state[i]]
        V_bleed[i] = 0.0
    else
        V_circ[i] = V[length(V)]
        V_bleed[i] = β[state[i] - length(V) + 1] * V[length(V)]
    end
end

# Compute the cut heating power
Q_cut = Q .- Q_aux[1:365*24]

# Figure
fig = Figure(; size = (17 * 96 / 2.54, 18 * 96 / 2.54))

ax = Axis(fig[1, 1], xlabel = L"$t$ (d)", ylabel = L"$Q$ (kW)")
lines!(ax, t[1:365*24] / (3600 * 24), Q / 1000, linewidth = 2, label = "Heat loads")
lines!(ax, t[1:365*24] / (3600 * 24), Q_aux[1:365*24] / 1000, linewidth = 2, label = "Aux. loads")
lines!(ax, t[1:365*24] / (3600 * 24), Q_cut[1:365*24] / 1000, linewidth = 2, label = "Cut loads")
axislegend(ax; position = :rc)

ax = Axis(fig[1, 2], xlabel = L"$t$ (y)", ylabel = L"$g$ (°Cm/W)", xscale = log10)
for gi in eachcol(g)
    lines!(ax, t / (3600 * 24 * 365), gi, linewidth = 2)
end

ax = Axis(fig[2, 1:2], xlabel = L"$t$ (y)", ylabel = L"$V$ (m³/s)")
lines!(ax, t / (3600 * 24 * 365), V_circ, linewidth = 2, label = "Circulating Flow Rate")
lines!(ax, t / (3600 * 24 * 365), V_bleed, linewidth = 2, label = "Bleed Flow Rate")
axislegend(ax; position = :rb)

ax = Axis(fig[3, 1:2], xlabel = L"$t$ (y)", ylabel = L"$T$ (°C)")
lines!(ax, t / (3600 * 24 * 365), T_in,  linewidth = 2, label = "Inlet Temp.")
lines!(ax, t / (3600 * 24 * 365), T_out, linewidth = 2, label = "Outlet Temp.")
hlines!(ax, Tlim, color = :black, linestyle = :dash, linewidth = 2)
axislegend(ax; position = :rb)

display(fig);
