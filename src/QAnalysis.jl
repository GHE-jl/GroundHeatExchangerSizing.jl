"""
Module that allows to construct a structure of ground thermal loads that can be used for further
    sizing or conception of ground source heat pump system.
"""

# Write your package code here.

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
        range = h_month[ii]+1:h_month[ii+1]
        Qm_ave[ii] = mean(Q[range])
        Qm_c_peak[ii] = minimum(Q[range])
        Qm_h_peak[ii] = maximum(Q[range])
    end

    # Correction
    Qm_c_peak[Qm_c_peak .> 0] .= 0.
    Qm_h_peak[Qm_h_peak .< 0] .= 0.

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

function building_to_ground_loads(Qb::Vector{T1}, COPh::T2, COPc::T2,
    pc_peakh=1, pc_peakc=1) where {T1<:Real, T2<:Real}
    """
    Converts building loads to ground loads, as a function of COP for both heating (COPh) and 
    cooling (COPc). Option inputs can be used to specify the percentage (pc) of the peak loads (for
    both heating and cooling) that has to be covered by the geothermal system. The default is 100%
    coverage.
    """
    #TODO Validate if it works

    # Cut the loads to the percentage of peak coverage.
    Qb[Qb .< pc_peakh*minimum(Qb)] .= pc_peakh*minimum(Qb)
    Qb[Qb .> pc_peakc*maximum(Qb)] .= pc_peakc*maximum(Qb)

    # Convert building loads (Qb) to ground loads (Qg)
    Qgh = Qb[Qb .<= 0] * (1 - (1/COPh))
    Qgc = Qb[Qb .> 0] * (1 + (1/COPc))
    Qg = Qgh + Qgc
    return Qg
end