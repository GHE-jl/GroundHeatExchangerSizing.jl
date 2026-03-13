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
    rmse = sqrt(sum(d .^ 2) / n)    # Root mean square error

    # Relative deviations; you can choose how to handle zeros in y
    rel = d ./ y                    # Relative deviation
    arel = abs.(rel)                # Absolute relative deviation
    mre  = sum(rel) / n             # Mean relative error
    mare = sum(arel) / n            # Mean absolute relative deviation

    return (mae = mae, mre = mre, mare = mare, rmse = rmse)
end