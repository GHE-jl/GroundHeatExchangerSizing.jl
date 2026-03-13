struct QLoads
    hr::Vector{Float64}         # Hourly annual heat load profile
    hₗ::Float64                 # Lower annual heat load value
    hₕ::Float64                 # Higher annual heat load value
    h::Float64                  # Absolute maximal heat load value
    m̄ₐ::Vector{Float64}         # Average monthly heat load values (12 values)
    m̄ₗ::Vector{Float64}         # Lower monthly loads (12 values)
    m̄ₕ::Vector{Float64}         # Higher monthly loads (12 values)
    mₗ::Float64                 # Average monthly heat load when Q.hₗ occurs
    mₕ::Float64                 # Average monthly heat load when Q.hₕ occurs
    y::Float64                  # Yearly average load
end

"""
    Q_analysis(Q)

Function that summarizes a thermal load vector to retrieve yearly, monthly and hourly peaks values
for both heating and cooling. These values can then be used to perform the sizing of a GHE with
different methods. The function returns a structure with a comprehensive summary of the heat load
profile.
# Arguments
    - Q: Heat load profile to analyse (8760x1) [W]
# Ouput
    - A structure QLoads based on the input vector of heat loads
"""
function Q_analysis(Q::AbstractVector{<:Real})
    # Preallocation
    Qm_ave = Vector(undef, 12)
    Qm_c_peak = Vector(undef, 12)
    Qm_h_peak = Vector(undef, 12)

    # Monthly heat load
    h_month = cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * 24)

    # Find values through months
    for ii in 1:12
        range = (h_month[ii] + 1):h_month[ii + 1]
        Qm_ave[ii] = sum(Q[range]) / length(Q[range])
        Qm_c_peak[ii] = minimum(Q[range])
        Qm_h_peak[ii] = maximum(Q[range])
    end

    # Correction
    Qm_c_peak[Qm_c_peak .> 0] .= 0.0
    Qm_h_peak[Qm_h_peak .< 0] .= 0.0

    return QLoads(
        Q,
        minimum(Q),
        maximum(Q),
        maximum(abs.(Q)),
        Qm_ave,
        Qm_c_peak,
        Qm_h_peak,
        Qm_ave[Qm_c_peak .== minimum(Q)][1],
        Qm_ave[Qm_h_peak .== maximum(Q)][1],
        sum(Q) / length(Q))
end

"""
    Q_hourly_to_monthly(Q)

Function that resample a hourly thermal load profile to monthly average and monthly peak loads
usually used in a monthly sizing equation. The monthly values have 36 outputs per year: 12 monthly 
average, 12 monthly cooling peak and 12 monthly heating loads.
# Argument
    - `Q`: Hourly thermal load profile (8760x1) [W]
# Output
    - Matrix of monthly thermal loads formed of [Qma, Qmc, Qmh] (12x3) [W]
"""
function Q_hourly_to_monthly(Q::AbstractVector{<:Real})
    # Monthly hours
    h_month = cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * 24)

    # Find values through months
    Qma, Qmc, Qmh = Vector(undef, 12), Vector(undef, 12), Vector(undef, 12)
    for i in 1:12
        range = (h_month[i] + 1):h_month[i + 1]
        Qma[i] = sum(Q[range]) / length(Q[range])
        Qmc[i] = minimum(Q[range])
        Qmh[i] = maximum(Q[range])
    end

    # Correction
    Qmc[Qmc .> 0] .= 0.0
    Qmh[Qmh .< 0] .= 0.0
    return [Qma Qmc Qmh]
end

"""
    Q_hourly_to_three_pulses(Q)

Function that resample a hourly thermal load profile to three pulses loads usually used in a
three pulses sizing equation. The three pulses loads are the yearly average `Qa`, the monthly average
when the peak load occurs `Qm` and the peak load `Qp`. Both `Qm` and `Qp` are retrieved for cooling
and heating.
# Argument
    - `Q`: Ground hourly thermal load profile (8760x1) [W]
# Output
    - Rows of three pulses loads (Qa; Qm; Qmp) and columns for cooling and heating (3x2) [W]
"""
function Q_hourly_to_three_pulses(Q::AbstractVector{<:Real})
    Qa = sum(Q) / length(Q)
    Qc = minimum(Q)
    Qh = maximum(Q)

    # Hourly to monthly
    Qm = Q_hourly_to_monthly(Q)
    Qmc = Qm[Qm[:, 2] .== Qc, 1]
    Qmh = Qm[Qm[:, 3] .== Qh, 1]

    return [Qa Qa; Qmc Qmh; Qc Qh]
end

"""
    Q_monthly_to_hourly(Q, tp)
    
Function that resample the monthly average and monthly peak thermal loads usually used in a monthly
sizing equation to hourly values for a year (used in hourly sizing equations). Monthly values have
36 inputs per year: 12 monthly average `Qma`, 12 monthly cooling peak `Qmc` and 12 monthly heating 
peak `Qmh` loads. This function converts these 36 loads to 8760 step loads considering peak time 
load at the end of each month.
# Arguments
    - `Q`: Matrix of monthly thermal loads formed of [Qma, Qmc, Qmh] (12x3) [W]
    - `tp`: Peak time (usually 4 to 6 hours) [s]
# Output
    - Cooling (first column) and heating (second column) hourly thermal loads (8760 x 2) [W]
"""
function Q_monthly_to_hourly(Q::AbstractMatrix{<:Real}, tp::Real)
    # Set all the time at which the load will change for the monthly method
    sec_month = repeat([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * 24 * 3600)
    sec_month_cumul = cumsum(sec_month)
    sec_month_peak = sec_month_cumul .- tp
    sec_id = sort([sec_month_cumul; sec_month_peak])

    # Inline the peak loads between the monthly average loads for cooling and heating
    Qvec_c, Qvec_h = zeros(24), zeros(24)   # Preallocation
    for i in 1:12                           # Cooling loads
        Qvec_c[2 * i - 1] = Q[i, 1]
        Qvec_c[2 * i] = Q[i, 2]
    end
    for i in 1:12                          # Heating loads
        Qvec_h[2 * i - 1] = Q[i, 1]
        Qvec_h[2 * i] = Q[i, 3]
    end

    # Set the correct heat load (average or peak) in a full length vector
    t = 3600:3600:(365 * 24 * 3600)
    Q = zeros(length(t), 2)

    for i in eachindex(sec_id)
        if i == 1
            range = t .<= sec_id[i]
        else
            range = (t .<= sec_id[i]) .& (t .> sec_id[i - 1])
        end
        Q[range, 1] .= Qvec_c[i]
        Q[range, 2] .= Qvec_h[i]
    end
    return Q
end

"""
    Q_monthly_to_three_pulses(Q)

Function that resample monthly thermal loads to three pulses, usually used in a three pulses sizing
equation. The three pulses loads are the yearly average `Qa`, the monthly average when the peak load
occurs `Qm` and the peak load `Qp`. Both `Qm` and `Qp` are retrieved for cooling and heating based
on monthly loads.
# Argument
    - `Q`: Matrix of monthly thermal loads formed of [Qma, Qmc, Qmh] (12x3) [W]
# Output
    - Rows of three pulses loads (Qa; Qm; Qmp) and columns for cooling and heating (3x2) [W]
"""
function Q_monthly_to_three_pulses(Q::AbstractMatrix{<:Real})
    Qa = sum(Q[:, 1]) / length(Q[:, 1])
    Qc = minimum(Q[:, 2])
    Qh = maximum(Q[:, 3])

    Qmc = Q[Q[:, 2] .== minimum(Q), 1]
    Qmh = Q[Q[:, 3] .== maximum(Q), 1]

    return [Qa Qa; Qmc Qmh; Qc Qh]
end