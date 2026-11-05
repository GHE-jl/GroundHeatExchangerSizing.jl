"""
    deviation_metrics(y, ŷ)

Function that computes deviation metrics between two vectors.
# Arguments
    - `y`: Observed signal.
    - `ŷ`: Predicted signal
# Outputs
    - `mae`: Mean absolute deviation
    - `mre`: Mean relative deviation (signed)
    - `mare`: Mean absolute relative deviation
    - `mse`: Mean squared deviation
    - `rmse`: Root mean squared deviation
"""
function deviation_metrics(y::AbstractVector, ŷ::AbstractVector)
    # Initial setup
    @assert length(y) == length(ŷ) "y and ŷ must have same length"
    n = length(y)

    # Difference between vectors
    d = ŷ .- y                      # Deviation between the vectors
    ae = abs.(d)                    # Absolute deviation

    # Basic deviation metrics
    mae  = sum(ae) / n              # Mean absolute error
    mse  = sum(d .^ 2) / n          # Mean squared error
    rmse = sqrt(mse)                # Root mean square error

    # Relative deviations; you can choose how to handle zeros in y
    rel = d ./ y                    # Relative deviation
    arel = abs.(rel)                # Absolute relative deviation
    mre  = sum(rel) / n             # Mean relative error
    mare = sum(arel) / n            # Mean absolute relative deviation

    return (mae = mae, mre = mre, mare = mare, mse = mse, rmse = rmse)
end

"""
    _Tlim_offset(Q, nb, V, ρf, cf, Tlim_at; loc = :mean)

Temperature offset [°C]  from the fluid temperature at the location `loc` to the fluid temperature
at `Tlim_at` (where the operating limits `Tlim` apply), for the ground load `Q` [W] (scalar or
vector) of `nb` boreholes with a flow rate `V` per borehole. The locations are seen from the ground,
not the heat pump. `:inlet` is the heat-pump leaving fluid temperature (LFT), `:outlet` the
heat-pump entering fluid temperature (EFT) and `:mean` the mean fluid temperature. With `Q > 0` for
heat injection (building cooling), the inlet is `Q/(2·nb·V·ρf·cf)` above the mean and the outlet the
same amount below it. `alternative_sizing` uses it from the mean fluid temperature (`loc = :mean`)
and `outlet_sizing` from the outlet temperature (`loc = :outlet`).
"""
function _Tlim_offset(Q, nb, V, ρf, cf, Tlim_at::Symbol; loc::Symbol = :mean)
    # Position of a location relative to the mean fluid temperature [loop ΔT]
    position(at) = at === :inlet ? 0.5 : at === :outlet ? -0.5 : at === :mean ? 0.0 :
        throw(ArgumentError("loc must be :inlet, :outlet or :mean, got :$at"))

    # Loop temperature difference [°C]
    ΔT = Q ./ (nb * V * ρf * cf)

    # Return the temperature offset from `loc` to `Tlim_at`
    return (position(Tlim_at) - position(loc)) .* ΔT
end

"""
    Q_COP(Qb, COP_heating, COP_cooling)

Convert building thermal load(s) to ground thermal load(s) using a heat-pump coefficient of
performance, so that a user given only a **building** load and average heating/cooling COPs can
produce the **ground** loads the sizing equations require. Inputs may be scalar or vector; the
GHE-jl sign convention is used:
- `Qb > 0` (building heating) → `Qg < 0` (heat extracted from the ground loop);
- `Qb < 0` (building cooling) → `Qg > 0` (heat rejected to the ground loop).

For heating, `Qg = -Qb·(1 - 1/COP_heating)`; for cooling, `Qg = -Qb·(1 + 1/COP_cooling)`.

This mirrors the `Q_COP` of `GroundSourceHeatPumpDesign.jl`; it is provided here as a convenience
access point. For load- or temperature-dependent COPs and capacity limiting, use that package.

# Arguments
    - `Qb`: Building thermal load, scalar or vector [W]
    - `COP_heating`: Heating COP, scalar (applied to every element) or per-element vector
    - `COP_cooling`: Cooling COP, scalar or per-element vector
# Output
    - `Qg`: Ground thermal load, same shape as `Qb` [W]
"""
function Q_COP(Qb::Real, COP_heating::Real, COP_cooling::Real)
    COP_heating > 0 || throw(ArgumentError("COP_heating must be strictly positive."))
    COP_cooling > 0 || throw(ArgumentError("COP_cooling must be strictly positive."))

    if Qb > 0
        return -Qb * (1 - 1 / COP_heating)
    elseif Qb < 0
        return -Qb * (1 + 1 / COP_cooling)
    end
    return zero(promote_type(Float64, typeof(Qb), typeof(COP_heating), typeof(COP_cooling)))
end

function Q_COP(Qb::AbstractVector{<:Real}, COP_heating::Real, COP_cooling::Real)
    T = promote_type(Float64, eltype(Qb), typeof(COP_heating), typeof(COP_cooling))
    Qg = similar(Qb, T)
    for k in eachindex(Qb)
        Qg[k] = Q_COP(Qb[k], COP_heating, COP_cooling)
    end
    return Qg
end

function Q_COP(Qb::AbstractVector{<:Real}, COP_heating::AbstractVector{<:Real},
    COP_cooling::AbstractVector{<:Real})
    length(Qb) == length(COP_heating) ||
        throw(DimensionMismatch("Qb and COP_heating must have the same length."))
    length(Qb) == length(COP_cooling) ||
        throw(DimensionMismatch("Qb and COP_cooling must have the same length."))
    T = promote_type(Float64, eltype(Qb), eltype(COP_heating), eltype(COP_cooling))
    Qg = similar(Qb, T)
    for k in eachindex(Qb, COP_heating, COP_cooling)
        Qg[k] = Q_COP(Qb[k], COP_heating[k], COP_cooling[k])
    end
    return Qg
end