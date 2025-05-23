"""
Module that allows to construct a structure of ground thermal loads that can be used for further
    sizing or conception of ground source heat pump system.
"""

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

function Q_analysis(Q)
    """
        Q_analysis(Q)

    Function that analyse a ground thermal load vector to retrieve yearly, monthly and hourly peaks 
    for both heating and cooling. These information are used primarily for ground heat exchanger 
    sizing. The convention used for this script is that building heating corresponds to ground
    cooling (so negative loads) and building cooling is ground heating (positive loads).
    Input:
        - Q: Heat load profile to analyse (8760x1) [W]
    Ouput:
        - A structure QLoads based on the input vector of heat loads
    """

    # Preallocation
    Qm_ave = Vector{Float64}(undef, 12)
    Qm_c_peak = Vector{Float64}(undef, 12)
    Qm_h_peak = Vector{Float64}(undef, 12)

    # Monthly heat load
    h_month = cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * 24)

    # Find values through months
    for ii in 1:12
        range = (h_month[ii] + 1):h_month[ii + 1]
        Qm_ave[ii] = mean(Q[range])
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
        mean(Q))
end

function Q_L3_resample(Q::Matrix{A}, t_peak::A, n_year) where A<:Real
    """
    Function that resample the monthly average and monthly peak loads usually used in a L3 sizing 
    equation and sample them to be used in an L4 method. More precisely, L3 methods needs 36 values 
    per year: 12 monthly average, 12 monthly cooling peak and 12 monthly heating loads. This 
    function converts these 36 loads to 8760*ny loads for cooling and heating considering a peak
    time load per month.
    Inputs:
        - Q: Matrix formed of [Q.m̄ₐ, Q.m̄ₗ, Q.m̄ₕ] from the QLoads structure (12x3) [W]
        - t_peak: Peak time (usually 4 to 6 hours) [s]
        - n_year: number of years used in simulation (usually 10) [-]
    Output:
        - Q: [8760*ny x 2] Array composed of:
            - First column: Cooling loads from monthly ones [W]
            - Second column: Heat loads from monthly ones [W]
    """

    # 1. Set all the time at which the load will change for the L3 method
    sec_month = repeat([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * 24 * 3600.0, Int(n_year))
    sec_month_cumul = cumsum(sec_month)
    sec_month_peak = sec_month_cumul .- (t_peak * 3600)
    sec_id = sort([sec_month_cumul; sec_month_peak])

    # 2. Inline the peak loads between the monthly average loads for cooling and heating
    Q_tmp = zeros(24)                   # Preallocation
    for ii in 1:12                      # Cooling loads
        Q_tmp[2 * ii - 1] = Q[ii, 1]
        Q_tmp[2 * ii] = Q[ii, 2]
    end
    Qvec_c = repeat(Q_tmp, Int(n_year))      # Repeat for the number of years

    for jj in 1:12                      # Heating loads
        Q_tmp[2 * jj - 1] = Q[jj, 1]
        Q_tmp[2 * jj] = Q[jj, 3]
    end
    Qvec_h = repeat(Q_tmp, Int(n_year))      # Repeat for the number of years

    # 3. Set the correct heat load (average or peak) in a full length vector
    t = range(1, 8760 * n_year, length = 8760 * Int(n_year)) * 3600
    Q = zeros(length(t), 2)

    for kk in eachindex(sec_id)
        if kk == 1
            range = t .<= sec_id[kk]
        else
            range = (t .<= sec_id[kk]) .& (t .> sec_id[kk - 1])
        end
        Q[range, 1] .= Qvec_c[kk]
        Q[range, 2] .= Qvec_h[kk]
    end

    return Q
end