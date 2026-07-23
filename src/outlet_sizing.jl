# Borehole-outlet transfer-function sizing (Dion & Pasquier 2025)
#
# The outlet temperature is the temporal superposition of the ground loads with a dimensionless
# borehole-outlet transfer function ḡ (`outlet_transfer_function` from GroundHeatExchanger.jl):
#
#     Tout(t) = Tg + Σᵢ (Qᵢ - Qᵢ₋₁)/(V·Cf) · ḡ(t - tᵢ₋₁, H)
#
# There is no explicit borehole resistance term, as Rb* is carried inside ḡ. Unlike the alternative
# ASHRAE equation (a fixed-point iteration), the borehole length here is found by optimisation. For
# each operating limit the length is the minimiser of |Tlim - extremum(Tout)| over H ∈ [50, 250] m
# — the H training range of `DeepANN`, the short-term ANN `outlet_transfer_function` defaults to
# (GroundHeatExchanger.jl).
#
# Ground loads Q are per convention negative for heat extraction (ground cooling, building heating)
# and positive for heat rejection (ground heating, building cooling).
#
# Reference:
#   - Dion, G., & and Pasquier, P. (2025). Ground heat exchanger sizing using borehole outlet 
#       transfer function. Science and Technology for the Built Environment, 31(10), 1–13. 
#       https://doi.org/10.1080/23744731.2025.2523200
#   - Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of artificial neural networks to
#       near-instant construction of short-term g-functions. Applied Thermal Engineering, 143,
#       910–921. https://doi.org/10.1016/j.applthermaleng.2018.07.137
#   - Pasquier, P., & Marcotte, D. (2020). Robust identification of volumetric heat capacity and
#       analysis of thermal response tests by Bayesian inference with correlated residuals. Applied
#       Energy, 261, 114394. https://doi.org/10.1016/j.apenergy.2019.114394

using Optim

# Bounded length-optimisation controls (not user choices).
const _H_LB = 50.0          # lower bound [m] (DeepANN H validity — the default short-term ANN)
const _H_UB = 250.0         # upper bound [m] (DeepANN H validity — the default short-term ANN)
const _H_TOL = 1.0e-3       # length convergence tolerance [m]
const _H_MAXIT = 1000       # optimiser iteration cap

"""
    _optimize_length(objective)

Minimise `objective(H)`: the temperature-limit residual `|Tlim - extremum(Tout(H))|` over the
borehole length `H ∈ [50, 250] m` (`DeepANN`'s H range). The residual is usually a smooth,
one-dimensional, unimodal (interior-minimum) function of `H`, so it is minimised with Optim.jl's
`Brent` method. Brent is a derivative-free bounded solver. It brackets the minimum inside
`[50, 250]` and shrinks the bracket until the length stops moving by more than `_H_TOL`. It needs
no initial guess, no gradient and no bound barrier, so it converges in a few tens of evaluations.
When a limit is never actually binding across the whole range (e.g. the low/heating limit for a
strongly cooling-dominated load), the residual instead decreases monotonically toward one of the
two edges, with its minimum sitting *at* that boundary rather than in the interior — a shape Brent
isn't built to find, since it only brackets interior dips. The two boundary values are therefore
checked explicitly against Brent's interior candidate, and whichever is smallest wins. Returns the
optimal length [m]. Internal helper shared by the three outlet levels.
"""
function _optimize_length(objective)
    res = optimize(objective, _H_LB, _H_UB, Brent(); abs_tol = _H_TOL, iterations = _H_MAXIT)
    H_interior, f_interior = Optim.minimizer(res), Optim.minimum(res)
    f_lo, f_hi = objective(_H_LB), objective(_H_UB)
    f_lo < f_interior && return _H_LB
    f_hi < f_interior && return _H_UB
    return H_interior
end

"""
    outlet_sizing_L2(Q3, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim;
        tp = 6.0, ny = 10.0)

Size a ground heat exchanger with the **three-pulse (L2)** borehole-outlet transfer-function
equation (Dion & Pasquier 2025): the outlet temperature after the yearly, monthly and peak
pulses is built from the transfer function evaluated at the three superposition times, and the
length is optimised against each operating limit.
# Arguments
    - `Q3`: Three-pulse ground loads (3 × 2) `[Qy Qy; Qm_c Qm_h; Qh_c Qh_h]` [W], column 1 cooling,
      column 2 heating — as returned by [`Q_hourly_to_three_pulses`](@ref)
    - `xy`: Borehole coordinates (nb × 2) [m]; `[0.0 0.0]` for a single borehole
    - `rb`: Borehole radius [m]
    - `D`: Buried depth [m]
    - `ks`, `Cs`: Ground thermal conductivity [W/mK] and volumetric heat capacity [J/m³K]
    - `s`: Shank spacing [m]
    - `ro`, `ri`: Pipe outer / inner radius [m]
    - `kg`, `Cg`: Grout thermal conductivity [W/mK] and volumetric heat capacity [J/m³K]
    - `kp`, `Cp`: Pipe thermal conductivity [W/mK] and volumetric heat capacity [J/m³K]
    - `kf`, `cf`, `ρf`, `μf`: Fluid thermal conductivity [W/mK], specific heat [J/kgK],
        density [kg/m³] and dynamic viscosity [kg/m/s]
    - `V`: Volumetric flow rate in one U-tube loop (per borehole) [m³/s]
    - `T0`: Undisturbed ground temperature [°C]
    - `Tlim`: Operating temperature limits `[low, high]` [°C]
# Keywords
    - `tp`: Peak pulse duration [h] (default 6)
    - `ny`: Design period [years] (default 10)
# Output
    - `H`: governing borehole length [m]
"""
function outlet_sizing_L2(Q3::AbstractMatrix{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, Cg::Real, kp::Real,
    Cp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, V::Real, T0::Real,
    Tlim::AbstractVector{<:Real}; tp::Real = 6.0, ny::Real = 10.0)
    nb = size(xy, 1)
    Cf = ρf * cf                                                 # fluid volumetric heat capacity
    ty = ny * 365 * 24 * 3600.0
    tm = round((((7 * 31) + (4 * 30) + 28) / 12) * 24 * 3600.0)
    th = tp * 3600.0
    tsup = [th, tm + th, ty + tm + th]                           # superposition times t1,t2,t3

    # Outlet temperature after the three pulses, optimised against each operating limit in turn.
    H_low = _optimize_length() do H
        Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
        ḡ = outlet_transfer_function(tsup, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V, s, Rbe,
            FLSModel(H, D, ks, Cs); xy = xy, model = DeepANN())
        Γh, Γm, Γy = ḡ[1], ḡ[2] - ḡ[1], ḡ[3] - ḡ[2]
        Tout = T0 + (Q3[1, 1] * Γy + Q3[2, 1] * Γm + Q3[3, 1] * Γh) / (nb * V * Cf)
        return abs(Tlim[1] - Tout)
    end
    H_high = _optimize_length() do H
        Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
        ḡ = outlet_transfer_function(tsup, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V, s, Rbe,
            FLSModel(H, D, ks, Cs); xy = xy, model = DeepANN())
        Γh, Γm, Γy = ḡ[1], ḡ[2] - ḡ[1], ḡ[3] - ḡ[2]
        Tout = T0 + (Q3[1, 2] * Γy + Q3[2, 2] * Γm + Q3[3, 2] * Γh) / (nb * V * Cf)
        return abs(Tlim[2] - Tout)
    end
    return max(H_low, H_high)
end

"""
    _outlet_convolution(Qc, Qh, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V,
        T0, Tlim, ny)

Shared convolution-based optimisation for the L3 and L4 outlet equations. `Qc`/`Qh` are the one-year
hourly ground load profiles [W] governing the low (cooling) and high (heating) limits, repeated over
`ny` years and superimposed with the outlet transfer function by `convolution` (Eqs. 10-12).
"""
function _outlet_convolution(Qc::AbstractVector{<:Real}, Qh::AbstractVector{<:Real},
    xy::AbstractMatrix{<:Real}, rb::Real, D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real,
    kg::Real, Cg::Real, kp::Real, Cp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, V::Real,
    T0::Real, Tlim::AbstractVector{<:Real}, ny::Real)
    nb = size(xy, 1)
    Cf = ρf * cf
    nyi = Int(ny)
    t = collect(3600.0:3600.0:(3600.0 * 8760 * nyi))
    Qc_full = repeat(collect(Float64, Qc), nyi)
    Qh_full = repeat(collect(Float64, Qh), nyi)

    # Hourly outlet temperature for a one-year profile repeated over the design period, optimised
    # against each operating limit in turn.
    H_low = _optimize_length() do H
        Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
        ḡ = outlet_transfer_function(t, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V, s, Rbe,
            FLSModel(H, D, ks, Cs); xy = xy, interp = true, model = DeepANN())
        Tout = T0 .+ convolution(Qc_full ./ (nb * V * Cf), ḡ)
        return abs(Tlim[1] - minimum(Tout))
    end
    H_high = _optimize_length() do H
        Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
        ḡ = outlet_transfer_function(t, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V, s, Rbe,
            FLSModel(H, D, ks, Cs); xy = xy, interp = true, model = DeepANN())
        Tout = T0 .+ convolution(Qh_full ./ (nb * V * Cf), ḡ)
        return abs(Tlim[2] - maximum(Tout))
    end
    return max(H_low, H_high)
end

"""
    outlet_sizing_L3(Qm, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim;
        tp = 6.0, ny = 10.0)

Size a ground heat exchanger with the **monthly (L3)** outlet transfer-function equation. The 36
monthly loads are expanded to hourly (monthly averages with the peak over the final `tp` hours of
each month, [`Q_monthly_to_hourly`](@ref)) and superimposed with the outlet transfer function.
# Arguments
    - `Qm`: Monthly ground loads (12 × 3) `[Qma Qmc Qmh]` [W] — as returned by
      [`Q_hourly_to_monthly`](@ref)
    - remaining arguments and keywords: as in [`outlet_sizing_L2`](@ref)
# Output
    - `H`: governing borehole length [m]
"""
function outlet_sizing_L3(Qm::AbstractMatrix{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, Cg::Real, kp::Real,
    Cp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, V::Real, T0::Real,
    Tlim::AbstractVector{<:Real}; tp::Real = 6.0, ny::Real = 10.0)
    Qhr = Q_monthly_to_hourly(Qm, tp * 3600.0)                   # 8760 × 2 [cooling, heating]
    return _outlet_convolution(Qhr[:, 1], Qhr[:, 2], xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp,
        kf, cf, ρf, μf, V, T0, Tlim, ny)
end

"""
    outlet_sizing_L4(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim;
        ny = 10.0)

Size a ground heat exchanger with the **hourly (L4)** outlet transfer-function equation,
superimposing the full 8760-hour ground load profile with the outlet transfer function by FFT
convolution (no load aggregation).
# Arguments
    - `Q`: Hourly ground load profile for one year (8760) [W]
    - remaining arguments: as in [`outlet_sizing_L2`](@ref) (no `tp`)
# Output
    - `H`: governing borehole length [m]
"""
function outlet_sizing_L4(Q::AbstractVector{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real, D::Real,
    ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, Cg::Real, kp::Real, Cp::Real,
    kf::Real, cf::Real, ρf::Real, μf::Real, V::Real, T0::Real, Tlim::AbstractVector{<:Real};
    ny::Real = 10.0)
    return _outlet_convolution(Q, Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf,
        V, T0, Tlim, ny)
end

"""
    outlet_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf, V, T0, Tlim;
        level = :L4, tp = 6.0, ny = 10.0)

Convenience entry point to the borehole-outlet transfer-function sizing taking the **hourly** ground
load profile `Q` (8760) [W] and resampling it internally for the requested `level`:
    - `:L2` (three pulses, [`Q_hourly_to_three_pulses`](@ref)),
    - `:L3` (monthly, [`Q_hourly_to_monthly`](@ref)),
    - `:L4` (hourly, default).
All other arguments and keywords are as in [`outlet_sizing_L2`](@ref).
# Output
    - `H`: governing borehole length [m]
"""
function outlet_sizing(Q::AbstractVector{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real, D::Real,
    ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, Cg::Real, kp::Real, Cp::Real,
    kf::Real, cf::Real, ρf::Real, μf::Real, V::Real, T0::Real, Tlim::AbstractVector{<:Real};
    level::Symbol = :L4, tp::Real = 6.0, ny::Real = 10.0)
    if level === :L2
        return outlet_sizing_L2(Q_hourly_to_three_pulses(Q), xy, rb, D, ks, Cs, s, ro, ri, kg, Cg,
            kp, Cp, kf, cf, ρf, μf, V, T0, Tlim; tp = tp, ny = ny)
    elseif level === :L3
        return outlet_sizing_L3(Q_hourly_to_monthly(Q), xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp,
            Cp, kf, cf, ρf, μf, V, T0, Tlim; tp = tp, ny = ny)
    elseif level === :L4
        return outlet_sizing_L4(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, Cg, kp, Cp, kf, cf, ρf, μf,
            V, T0, Tlim; ny = ny)
    else
        throw(ArgumentError("level must be :L2, :L3 or :L4, got :$level"))
    end
end
