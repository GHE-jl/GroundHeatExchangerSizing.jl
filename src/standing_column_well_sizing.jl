using Optim
includet("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

"""
    scw_sizing_three_pulses(Q, Tlim, tp, k, kp, rb, ro, ri, V, β, T₀)

Optimizes standing column well (SCW) length using a n-year monthly ground load profile so that
operating temperature remains between operational limits. This script assumes that the ground has
uniform thermal and hydraulic conductivity. The monthly version interpolates the loads to create
hourly profiles for convolution, with the peak loads applied for tp hours. The function is defined
for either a single or multiple flow rates.
# Arguments
    - `Q`: Rows of thermal loads (`Qa`; `Qm`; `Qp`) and columns for cooling and heating (3x2) [W]
        - `Qa`: Yearly average, `Qm`: Monthly average when peak load occurs, `Qp`: Peak load
    - `Tlim`: Lower and higher fluid operating temperature limit (2x1) [°C]
    - `tp`: Time of the peak heat pulse (typically either 4 or 6 hours) [h]
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
    - Dion, G., & Pasquier, P. (2025). Adaptation of sizing equations for standing column wells. 
        Proceedings of the European Geothermal Congress 2025, 8. 
        https://europeangeothermalcongress.eu/wp-content/uploads/2025/11/Dion-et-Philippe.pdf
"""
function scw_sizing_three_pulses(Q, Tlim, tp, ks, kp, Cs, K, rb, ro, ri, V::Real, β::Real, T₀)
    # 1. Time superposition and constant definition
    ty = 10*365*24*3600                             # Seconds in 10 years
    tm = round((((7*31)+(4*30)+28)/12)*24*3600)     # Average seconds in 1 month
    t = [tp, tp + tm, tp + tm + ty]                 # Time vector for the three pulses (s)
    Cf = water_cp(T₀) * water_ρ(T₀)                 # (GHEModels.jl) Cf at T₀
    
    # 2. Define Objective Function
    function obj_fun(Hᵢ, limit_T, Q)
        # 2.1 Compute the SCW transfer function using the βILS model in GHEModels.jl
        g = βils_outlet(t, [Hᵢ ks Cs], kp, [Hᵢ K], rb, ro, ri, Hᵢ, Hᵢ-10.0, V, β, T₀)
        
        # 2.2 Compute T_out using convolution
        T_out = T₀ + (Q[1] * (g[3] - g[2]) + Q[2] * (g[2] - g[1]) + Q[3] * g[1] ) / (V * Cf)
        T_in = T_out .+ Q / (V * Cf)                # Inlet temperature from energy balance
        println("T_out: ", T_out, " °C, T_in: ", T_in, " °C")

        # 2.3 Evaluate objective based on Tlim type
        return abs(T_in - limit_T)
    end

    # 3. Run Optimization
    H_min, H_max = 100.0, 500.0
    H = zeros(2)
    out = Vector{Any}()
    for i in 1:2
        res = optimize(H -> obj_fun(H, Tlim[i], Q[:, i]),
            H_min, H_max, Brent(), iterations = 100, rel_tol = 1e-3, abs_tol = 1e-3)
        H[i] = Optim.minimizer(res)
        push!(out, res)
        
        println("Optimization $(i) (Tlim=$(Tlim[i])): H = $(round(H[i], digits=3)) m")
    end
    return maximum(H), out
end

"""
    scw_sizing_monthly(Q, Tlim, ny, tp, k, kp, rb, ro, ri, V, β, T₀)

Optimizes standing column well (SCW) length using a n-year monthly ground load profile so that
operating temperature remains between operational limits. This script assumes that the ground has
uniform thermal and hydraulic conductivity. The monthly version interpolates the loads to create
hourly profiles for convolution, with the peak loads applied for tp hours. The function is defined
for either a single or multiple flow rates.
# Arguments
    - `Q`: Matrix of monthly thermal loads formed of [Qma, Qmc, Qmh] (12x3) [W]
        - `Qma`: Monthly average, `Qmc`: Monthly cooling peak, `Qmh`: Monthly heating peak
    - `Tlim`: Lower and higher fluid operating temperature limit (2x1) [°C]
    - `ny`: Number of years to cover (typically 10) [-]
    - `tp`: Time of the peak heat pulse (typically either 4 or 6 hours) [h]
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
    - Dion, G., & Pasquier, P. (2025). Adaptation of sizing equations for standing column wells. 
        Proceedings of the European Geothermal Congress 2025, 8. 
        https://europeangeothermalcongress.eu/wp-content/uploads/2025/11/Dion-et-Philippe.pdf
"""
function scw_sizing_monthly(Q_, Tlim, ny, tp, ks, kp, Cs, K, rb, ro, ri, V::Real, β::Real, T₀)
    # 1. Thermal loads setup
    t = 3600.0:3600.0:(ny * 365 * 24 * 3600.0)
    t_ = t[set_nodes(length(t), 100)]               # From GHEModels.jl (Utils) reduced time vector
    Q = repeat(Q_monthly_to_hourly(Q_, tp), ny)     # Repeat the yearly profile for n years
    dQ = hcat([impulse_func(col) for col in eachcol(Q)]...) # Compute the impulse thermal loads
    Cf = water_cp(T₀) * water_ρ(T₀)                 # (GHEModels.jl) Cf at T₀
    
    # 2. Define Objective Function
    function obj_fun(Hᵢ, limit_T, Q, dQ)
        # 2.1 Compute the SCW transfer function using the βILS model in GHEModels.jl
        g_ = βils_outlet(t_, [Hᵢ ks Cs], kp, [Hᵢ K], rb, ro, ri, Hᵢ, Hᵢ-10.0, V, β, T₀)
        g = pchip_interpolation(t_, g_, t)          # From GHEModels.jl (Utils)
        
        # 2.2 Compute T_out using convolution
        T_out = T₀ .+ convolution(dQ / Hᵢ, g)
        T_in = T_out .+ Q / (V * Cf)                # Inlet temperature from energy balance

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
    out = Vector{Any}()
    for i in 1:2
        res = optimize(H -> obj_fun(H, Tlim[i], Q[:, i], dQ[:, i]),
            H_min, H_max, Brent(), iterations = 100, rel_tol = 1e-3, abs_tol = 1e-3)
        H[i] = Optim.minimizer(res)
        push!(out, res)
        
        println("Optimization $(i) (Tlim=$(Tlim[i])): H = $(round(H[i], digits=3)) m")
    end
    return maximum(H), out
end

"""
    scw_sizing_hourly(Q, Tlim, ny, k, kp, rb, ro, ri, V, β, T₀)

Optimizes standing column well (SCW) length using a n-year hourly ground load profile so that
operating temperature remains between operational limits. This script assumes that the ground has
uniform thermal and hydraulic conductivity. The function is defined for either a single or multiple
flow rates.
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
    - `V`: Fluid circulation flow rate (1x1 or nₛx1) [m³/s]
    - `β`: Fluid bleed rate ratio (between 0.01 and 1) (1x1 or nₛx1) [-]
    - `T₀`: Initial ground temperature [°C] (used for water properties)
# Outputs
    - `maximum(H)`: Maximal sizing length to cover the thermal loads [m]
    - `out`: Summary of the optimization outputs for both temperature limits [-]
# Reference
    - Dion, G., & Pasquier, P. (2025). Adaptation of sizing equations for standing column wells. 
        Proceedings of the European Geothermal Congress 2025, 8. 
        https://europeangeothermalcongress.eu/wp-content/uploads/2025/11/Dion-et-Philippe.pdf
"""
function scw_sizing_hourly(Q_, Tlim, ny, ks, kp, Cs, K, rb, ro, ri, V::Real, β::Real, T₀)
    # 1. Thermal loads setup
    t = 3600.0:3600.0:(ny * 365 * 24 * 3600.0)
    t_ = t[set_nodes(length(t), 100)]               # From GHEModels.jl (Utils) reduced time vector
    Q = repeat(Q_, ny)                              # Repeat the yearly profile for n years
    dQ = impulse_func(Q)                            # Compute the impulse thermal loads
    Cf = water_cp(T₀) * water_ρ(T₀)                 # (GHEModels.jl) Cf at T₀
    
    # 2. Define Objective Function
    function obj_fun(Hᵢ, limit_T)
        # 2.1 Compute the SCW transfer function using the βILS model in GHEModels.jl
        g_ = βils_outlet(t_, [Hᵢ ks Cs], kp, [Hᵢ K], rb, ro, ri, Hᵢ, Hᵢ-10.0, V, β, T₀)
        g = pchip_interpolation(t_, g_, t)          # From GHEModels.jl (Utils)
        
        # 2.2 Compute T_out using convolution
        T_out = T₀ .+ convolution(dQ / Hᵢ, g)
        T_in = T_out .+ Q / (V * Cf)                # Inlet temperature from energy balance

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
    out = Vector{Any}()
    for i in 1:2
        res = optimize(H -> obj_fun(H, Tlim[i]),
            H_min, H_max, Brent(), iterations = 100, rel_tol = 1e-3, abs_tol = 1e-3)
        H[i] = Optim.minimizer(res)
        push!(out, res)
        
        println("Optimization $(i) (Tlim=$(Tlim[i])): H = $(round(H[i], digits=3)) m")
    end
    return maximum(H), out
end

function scw_sizing_hourly(Q_, Tlim, ny, ks, kp, Cs, K, rb, ro, ri, V::AbstractVector,
    β::AbstractVector, T₀)
    # 1. Time setup
    t = 3600.0:3600.0:(ny * 365 * 24 * 3600.0)
    t_ = t[set_nodes(length(t), 100)]               # From GHEModels.jl (Utils) reduced time vector
    Q = repeat(Q_, ny)                              # Repeat the yearly profile for n years
    global g = zeros(length(t), length(V) + length(β) - 1)  # Preallocate for transfer functions
    
    # 2. Define Objective Function
    function obj_fun(Hᵢ, limit_T)
        # 2.1 Compute the SCW transfer function using the βILS model in GHEModels.jl
        for i in eachindex(V)
            if i != length(V)
                tmp = βils_outlet(t_, [Hᵢ ks Cs], kp, [Hᵢ K], rb, ro, ri, Hᵢ, Hᵢ-10.0, V[i], β[1],
                    T₀)
                g[:, i] = pchip_interpolation(t_, tmp, t) # From GHEModels.jl (Utils)
            else
                for j in eachindex(β)
                    tmp = βils_outlet(t_, [Hᵢ ks Cs], kp, [Hᵢ K], rb, ro, ri, Hᵢ, Hᵢ-10.0, V[i],
                        β[j], T₀)
                    g[:, i + j - 1] = pchip_interpolation(t_, tmp, t) # From GHEModels.jl (Utils)
                end
            end
        end

        # 2.2 Optimize the SCW operation
        state = ones(Int, length(t))                # Initial state (flow rate index)
        T_out = similar(t)
        T_in  = similar(t)

        e_aux, state, _, _, _ = states_optim!(T_in, T_out, state, Q/Hᵢ, g, V, T₀;
            mode = limit_T < 17.5 ? :heating : :cooling, ΔT = 0.3, maxiter = 25)

        # 2.3 Added constraint to the objective function
        state_penalty = sum(state .- 1)             # Penalize higher flow rates
        println("Auxiliary Energy: ", e_aux, " J, State Penalty: ", state_penalty)

        return e_aux + state_penalty + Hᵢ
    end

    # 3. Run Optimization
    H_min, H_max = 100.0, 500.0
    H = zeros(2)
    out = Vector{Any}()
    for i in 1:2
        res = optimize(H -> obj_fun(H, Tlim[i]),
            H_min, H_max, Brent(), iterations=100, rel_tol = 1e-3, abs_tol = 1e-3)
        H[i] = Optim.minimizer(res)
        push!(out, res)
        
        println("Optimization $(i) (Tlim=$(Tlim[i])): H = $(round(H[i], digits=3)) m")
    end
    return maximum(H), out
end

"""
    states_optim!(T_in, T_out, state, Q, g, V, T₀; mode = :heating, ΔT = 0.3, maxiter = 25)

Optimization of the flow rate state at each time step based on temperature thresholds for heating 
and cooling mode.
# Arguments
    - `T_in`: Vector to store the inlet temperature at each time step (modified in-place) [°C]
    - `T_out`: Vector to store the outlet temperature at each time step (modified in-place) [°C]
    - `state`: Vector of flow rate state indices (modified in-place) [-]
    - `Q`: Thermal load vector (ny*8760x1) [W/m]
    - `g`: Matrix of SCW transfer functions for each flow rate and bleed ratio combination [°Cm/W]
    - `V`: Vector of flow rates corresponding to the states [m³/s]
    - `T₀`: Initial ground temperature [°C] (used for water properties)
    - `mode`: Symbol for :heating or :cooling linked to the logic to use.
    - `ΔT` (default 0.3): Deadband temperature to prevent rapid state switching [°C]
    - `maxiter` (default 25): Maximum number of iterations to perform [-]
# Outputs
    - `auxiliary_energy`: Total energy added or removed by the auxiliary system to maintain 
        temperature limits [J]
    - `state`: Final vector of flow rate state indices after optimization [-]
"""
function states_optim!(T_in, T_out, state, Q, g, V, T₀; mode = :heating, ΔT = 0.3, maxiter = 25)

    # Set the thresholds (in °C) for state transitions based on mode
    if mode == :heating
        Tth = [10.0, 7.0, 5.0, 3.0, 1.0]
        T_aux = Tth[end]
    else
        Tth = [25.0, 28.0, 30.0, 32.0, 34.0]
        T_aux = Tth[end]
    end

    nt = length(T_in)
    e_aux = 0.0
    iter = 0
    changed = true

    while changed && iter < maxiter
        iter += 1
        old_state = copy(state)

        # Compute non-stationary convolution
        T_out .= T₀ .+ convolution_ns(Q, g, state)

        # Update state at each time step
        for k in 1:nt
            Cf = water_cp(T_out[k]) * water_ρ(T_out[k])
            # Compute T_in
            V_ = V[min(state[k], length(V))] # Ensure state index does not exceed max flow rates
            T_in[k] = T_out[k] - Q[k] / (V_ * Cf)

            # State selection based on mode, thresholds and temperature value
            new_state = select_state(T_in[k], state[k], Tth, ΔT, mode)

            # Auxiliary shaving
            if mode == :heating && T_in[k] ≤ T_aux
                deficit = (T_aux - T_in[k]) * V_ * Cf
                e_aux += deficit
                Q[k] -= deficit
                new_state = length(Tth)
            elseif  mode == :cooling && T_in[k] ≥ T_aux
                excess = (T_in[k] - T_aux) * V_ * Cf
                e_aux += excess
                Q[k] -= excess
                new_state = length(Tth)
            end

            # Smoothing of state transition: forbid sudden recovery
            if k > 1
                if mode == :heating
                    new_state = max(new_state, state[k-1])
                else
                    new_state = max(new_state, state[k-1])
                end
            end
            state[k] = new_state
        end
        # Keep going until no state change anymore or max iterations reached
        changed = any(state .!= old_state)
    end
    return e_aux, state, T_in, T_out, Q
end

"""
    select_state(T, prev_state, Tth, ΔT, mode)

Logic for state transition based on temperature thresholds for heating and cooling modes. The state
changes by no more than one increment by computation.
# Arguments
    - `T`: Temperature value [°C]
    - `prev_state`: State before checking for changes [-]
    - `Tth`: Vector of thresholds according to the `mode` (nₛ x 1) [°C]
        - mode = :heating, thresholds in descending ordre (e.g., [15.0, 10.0, 5.0])
        - mode = :cooling, thresholds in ascending order (e.g., [25.0, 30.0, 35.0])
    - `ΔT`: Deadband temperature to prevent rapid state switching [°C]
    - `mode`: Symbol for :heating or :cooling linked to the logic to use.
# Output
    - The new state value [-]
# Example
    select_state(12.5, 1, [15.0, 10.0, 5.0], 0.0, :heating)  # Output: new_state = 2,  <15.0 °C
    select_state(16.0, 3, [15.0, 10.0, 5.0], 0.0, :heating)  # Output: new_state = 2, Max change: 1
    select_state(5.1, 3, [15.0, 10.0, 5.0], 3.0, :heating)  # Output: new_state = 3, Because of ΔT
"""
function select_state(T::Float64, prev_state::Int, Tth::Vector{Float64}, ΔT::Float64, mode::Symbol)
    # Initialization
    n = length(Tth)         # Maximum number of state
    new_state = prev_state  # Consider no change in state

    # Building heating (ground cooling)
    if mode == :heating
        # Worsening condition (colder ground)
        if prev_state < n && T ≤ Tth[prev_state] - ΔT
            new_state = prev_state + 1
        end
        # Improving condition (warmer ground)
        if prev_state > 1 && T ≥ Tth[prev_state] + ΔT
            new_state = prev_state - 1
        end
    # Building cooling (ground heating)
    else
        # Worsening condition (warmer ground)
        if prev_state < n && T ≥ Tth[prev_state] + ΔT
            new_state = prev_state + 1
        end
        # Improving condition (cooler ground)
        if prev_state > 1 && T ≤ Tth[prev_state] - ΔT
            new_state = prev_state - 1
        end
    end
    return new_state
end