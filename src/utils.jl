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