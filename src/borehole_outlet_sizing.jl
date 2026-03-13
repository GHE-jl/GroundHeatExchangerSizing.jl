using Optim
include("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

"""
    outlet_sizing_hourly(Q, Tlim, ny, k, kp, rb, ro, ri, V, β, T₀)

Optimizes standing column well (SCW) length using a n-year hourly ground load profile so that
operating temperature remains between operational limits. This script assumes that the ground has
uniform thermal and hydraulic conductivity.
# Arguments
    - `Q`: Ground hourly thermal load profile (8760x1) [W]
    - `Tlim`: Lower and higher fluid operating temperature limit (2x1) [°C]
    - `ny`: Number of years to cover (typically 10) [-]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `K`: Groundwater hydraulic conductivity [m/s] 
    - `rb`: Borehole radius [m]
    - `ro`: Outer pipe radius [m]
    - `ri`: Inner pipe radius [m]
    - `V`: Fluid circulation flow rate [m³/s]
    - `β`: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - `T₀`: Initial ground temperature [°C] (used for water properties)
# Outputs
    - `maximum(H)`: Maximal sizing length to cover the thermal loads [m]
    - `out`: Summary of the optimization outputs for both temperature limits [-]
# Reference
    - ...
"""
function outlet_sizing_hourly(Q, Tlim, ny, ks, kp, Cs, K, rb, ro, ri, V, β, T₀)
    # 1. Thermal loads setup
    dQ = diff([0; repeat(Q, ny)])
    t = 3600.0:3600.0:(ny * 365 * 24 * 3600.0)
    t_ = t[set_nodes(length(t), 100) ]              # From GHEModels.jl (Utils) reduced time vector
    
    # 2. Define Objective Function
    function obj_fun(Hᵢ, limit_T, t, t_, ks, kp, Cs, K, rb, ro, ri, V, β, T₀, dQ)
        # 2.1 Compute the SCW transfer function using the βILS model in GHEModels.jl
        # g_ = Analytical model here!
        g = pchip_interpolation(t_, g_, t)          # From GHEModels.jl (Utils)
        
        # 2.2 Compute T_out using convolution
        T_out = T₀ .+ convolution(dQ / Hᵢ, g)
        T_in = T_out .- dQ / (V * Cs)               # Inlet temperature from energy balance

        # 2.3 Evaluate objective based on Tlim type
        if limit_T < 10
            # Ground cooling (Heating mode) -> minimize difference at minimum temp
            return abs(minimum(T_in) - limit_T)
        else
            # Ground heating (Cooling mode) -> minimize difference at maximum temp
            return abs(maximum(T_in) - limit_T)
        end
    end

    # 3. Run Optimization
    H_min, H_max = 100.0, 500.0
    H = zeros(2)
    # results = zeros(2)
    out = []
    for i in 1:2
        res = optimize(H -> obj_fun(H, Tlim[i], t, t_, ks, kp, Cs, K, rb, ro, ri, V, β, T₀, dQ),
            H_min, H_max, method = Brent(),
            iterations=100, rel_tol = 1e-3, abs_tol = 1e-3)
        H[i] = Optim.minimizer(res)
        push!(out, res)
        
        println("Optimization $(i) (Tlim=$(Tlim[i])): H = $(round(H[i], digits=3)) m")
    end
    return maximum(H), out
end