# Alternative ASHRAE sizing equation (g-function based)
#
# The "alternative" ASHRAE sizing equation of Ahmadfard & Bernier (2018, 2019) sizes a vertical
# ground heat exchanger by expressing the mean fluid temperature as the temporal superposition of
# the ground thermal loads with a finite-line-source (FLS) g-function evaluated for the actual
# borefield.
#
# Ground loads Q are per convention negative for heat extraction (ground cooling, building heating)
# and positive for heat rejection (ground heating, building cooling).
#
# Reference:
#   - Ahmadfard, M., & Bernier, M. (2018). Modifications to ASHRAE’s sizing method for vertical 
#       ground heat exchangers. Science and Technology for the Built Environment, 24(7), 803–817. 
#       https://doi.org/10.1080/23744731.2018.1423816
#   - Ahmadfard, M., & Bernier, M. (2019). A review of vertical ground heat exchanger sizing tools
#       including an inter-model comparison. Renewable and Sustainable Energy Reviews, 110, 247–265.
#       https://doi.org/10.1016/j.rser.2019.04.045

# Fixed-point iteration controls.
const _ALT_TOL = 0.01      # borehole-length convergence tolerance [m]
const _ALT_MAXIT = 20      # maximum number of iterations

"""
    alternative_sizing_L2(Q3, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim;
        tp = 6.0, ny = 10.0, Tlim_at = :inlet)

Size a ground heat exchanger with the **three-pulse (L2)** alternative ASHRAE equation. The load is
reduced to three ground pulses — yearly average, monthly average during the peak month, and the
peak, applied over `ny` years, one month, and `tp` hours. The three effective ground thermal
resistances come from the field FLS g-function evaluated at those superposition times, and the
borehole length is found by the fixed-point iteration of Ahmadfard & Bernier (2019).
# Arguments
    - `Q3`: Three-pulse ground loads (3 × 2) `[Qy Qy; Qm_c Qm_h; Qh_c Qh_h]` [W], column 1 cooling
        (drives the low limit), column 2 heating — as returned by [`Q_hourly_to_three_pulses`](@ref)
    - `xy`: Borehole coordinates (nb × 2) [m]; `[0.0 0.0]` for a single borehole
    - `rb`: Borehole radius [m]
    - `D`: Buried depth [m]
    - `ks`, `Cs`: Ground thermal conductivity [W/mK] and volumetric heat capacity [J/m³K]
    - `s`: Shank spacing [m]
    - `ro`, `ri`: Pipe outer / inner radius [m]
    - `kg`, `kp`: Grout and pipe thermal conductivity [W/mK]
    - `kf`, `cf`, `ρf`, `μf`: Fluid thermal conductivity [W/mK], specific heat [J/kgK],
        density [kg/m³] and dynamic viscosity [kg/m/s]
    - `V`: Volumetric flow rate in one U-tube loop (per borehole) [m³/s]
    - `T0`: Undisturbed ground temperature [°C]
    - `Tlim`: Operating temperature limits `[low, high]` [°C], applied at `Tlim_at`
# Keywords
    - `tp`: Peak pulse duration [h] (default 6)
    - `ny`: Design period [years] (default 10)
    - `Tlim_at`: Where `Tlim` applies, seen from the ground perspective: `:inlet` (ground
        inlet = heat-pump LFT, default, the most conservative), `:outlet` (ground outlet = heat-pump
        EFT) or `:mean` (mean fluid temperature)
# Output
    - `H`: governing borehole length [m]
"""
function alternative_sizing_L2(Q3::AbstractMatrix{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real,
    cf::Real, ρf::Real, μf::Real, V::Real, T0::Real, Tlim::AbstractVector{<:Real};
    tp::Real = 6.0, ny::Real = 10.0, Tlim_at::Symbol = :inlet)
    # Basic inputs
    _check_flow_rate(V)
    nb = size(xy, 1)
    offset = _Tlim_offset(Q3[3, :], nb, V, ρf, cf, Tlim_at; loc = :mean) # peak offset [°C]
    ty = ny * 365 * 24 * 3600.0                                  # design period [s]
    tm = round((((7 * 31) + (4 * 30) + 28) / 12) * 24 * 3600.0)  # average month [s]
    th = tp * 3600.0                                             # peak duration [s]
    tsup = [th, tm + th, ty + tm + th]                           # superposition times t1,t2,t3

    Hi = [0.0, 150.0]
    while abs(Hi[end] - Hi[end - 1]) > _ALT_TOL && length(Hi) < _ALT_MAXIT
        # Initiation of the iterative sizing loop
        H = Hi[end]

        # Per-borehole field g-function [°Cm/W] and field ground resistance.
        g = ground_response(tsup, rb, xy, FLSModel(H, D, ks, Cs); bc = :II, solver = :successive)
        Rgh, Rgm, Rgy = g[1], g[2] - g[1], g[3] - g[2]           # effective ground resistances
        Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

        # Length required to meet the low and high temperature limits
        L_low  = (Q3[1, 1] * Rgy + Q3[2, 1] * Rgm + Q3[3, 1] * (Rgh + Rbe)) /
            (Tlim[1] - T0 - offset[1])
        L_high = (Q3[1, 2] * Rgy + Q3[2, 2] * Rgm + Q3[3, 2] * (Rgh + Rbe)) /
            (Tlim[2] - T0 - offset[2])
        push!(Hi, max(L_low, L_high) / nb)
    end
    return Hi[end]
end

"""
    _alternative_convolution(Qc, Qh, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0,
        Tlim, ny, Tlim_at)

Shared convolution-based fixed-point sizing for the L3 and L4 alternative equations. `Qc`/`Qh` are
the one-year hourly ground load profiles [W] governing the low (cooling) and high (heating) limits;
they are repeated over `ny` years and superimposed with the field FLS g-function by `convolution`.
"""
function _alternative_convolution(Qc::AbstractVector{<:Real}, Qh::AbstractVector{<:Real},
    xy::AbstractMatrix{<:Real}, rb::Real, D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real,
    kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, V::Real, T0::Real,
    Tlim::AbstractVector{<:Real}, ny::Real, Tlim_at::Symbol)

    # Basic inputs
    _check_flow_rate(V)
    nb = size(xy, 1)
    nyi = Int(ny)
    t = collect(3600.0:3600.0:(3600.0 * 8760 * nyi))             # hourly time grid over ny years
    Qc_full = repeat(collect(Float64, Qc), nyi)
    Qh_full = repeat(collect(Float64, Qh), nyi)
    offsetc = _Tlim_offset(Qc_full, nb, V, ρf, cf, Tlim_at; loc = :mean)      # hourly offsets [°C]
    offseth = _Tlim_offset(Qh_full, nb, V, ρf, cf, Tlim_at; loc = :mean)

    Hi = [0.0, 150.0]
    while abs(Hi[end] - Hi[end - 1]) > _ALT_TOL && length(Hi) < _ALT_MAXIT

        # Initiation of the current iteration
        H = Hi[end]

        # Per-borehole field g-function [°Cm/W] → field effective ground resistance.
        g = ground_response(t, rb, xy, FLSModel(H, D, ks, Cs); bc = :II, solver = :successive)
        Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

        # Required borehole lengths for low and high limits based on convolution with the g-function
        L_low  = maximum((convolution(Qc_full, g) .+ Qc_full .* Rbe) ./ (Tlim[1] - T0 .- offsetc))
        L_high = maximum((convolution(Qh_full, g) .+ Qh_full .* Rbe) ./ (Tlim[2] - T0 .- offseth))
        push!(Hi, max(L_low, L_high) / nb)
    end
    return Hi[end]
end

"""
    alternative_sizing_L3(Qm, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim;
                          tp = 6.0, ny = 10.0, Tlim_at = :inlet)

Size a ground heat exchanger with the **monthly (L3)** alternative ASHRAE equation. The 36 monthly
loads are expanded to an hourly profile (monthly averages with the peak over the final `tp` hours of
each month, [`Q_monthly_to_hourly`](@ref)) and superimposed with the field FLS g-function.
# Arguments
    - `Qm`: Monthly ground loads (12 × 3) `[Qma Qmc Qmh]` [W] — as returned by
      [`Q_hourly_to_monthly`](@ref)
    - remaining arguments and keywords: as in [`alternative_sizing_L2`](@ref)
# Output
    - `H`: governing borehole length [m]
"""
function alternative_sizing_L3(Qm::AbstractMatrix{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real,
    cf::Real, ρf::Real, μf::Real, V::Real, T0::Real, Tlim::AbstractVector{<:Real};
    tp::Real = 6.0, ny::Real = 10.0, Tlim_at::Symbol = :inlet)

    # Convert monthly loads to hourly profile for convolution-based sizing
    Qhr = Q_monthly_to_hourly(Qm, tp * 3600.0)                   # 8760 × 2 [cooling, heating]

    # Perform convolution-based sizing using the hourly load profile
    return _alternative_convolution(Qhr[:, 1], Qhr[:, 2], xy, rb, D, ks, Cs, s, ro, ri, kg, kp,
        kf, cf, ρf, μf, V, T0, Tlim, ny, Tlim_at)
end

"""
    alternative_sizing_L4(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim;
                          ny = 10.0, Tlim_at = :inlet)

Size a ground heat exchanger with the **hourly (L4)** alternative ASHRAE equation, superimposing the
full 8760-hour ground load profile with the field FLS g-function by FFT convolution.
# Arguments
    - `Q`: Hourly ground load profile for one year (8760) [W]
    - remaining arguments and keywords: as in [`alternative_sizing_L2`](@ref) (no `tp`: hourly
      loads are used directly)
# Output
    - `H`: governing borehole length [m]
"""
function alternative_sizing_L4(Q::AbstractVector{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real,
    cf::Real, ρf::Real, μf::Real, V::Real, T0::Real, Tlim::AbstractVector{<:Real};
    ny::Real = 10.0, Tlim_at::Symbol = :inlet)

    # Perform convolution-based sizing using the hourly load profile for both cooling and heating
    return _alternative_convolution(Q, Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf,
        V, T0, Tlim, ny, Tlim_at)
end

"""
    alternative_sizing(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0, Tlim;
        level = :L4, tp = 6.0, ny = 10.0, Tlim_at = :inlet)

Convenience entry point to the alternative ASHRAE sizing equation taking the **hourly** ground load
profile `Q` (8760) [W] and resampling it internally for the requested `level`:
    - `:L2` (three pulses, [`Q_hourly_to_three_pulses`](@ref)),
    - `:L3` (monthly, [`Q_hourly_to_monthly`](@ref))
    - `:L4` (hourly, default).
All other arguments and keywords are as in [`alternative_sizing_L2`](@ref).
# Output
    - `H`: governing borehole length [m]
"""
function alternative_sizing(Q::AbstractVector{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real,
    cf::Real, ρf::Real, μf::Real, V::Real, T0::Real, Tlim::AbstractVector{<:Real};
    level::Symbol = :L4, tp::Real = 6.0, ny::Real = 10.0, Tlim_at::Symbol = :inlet)

    # Select the appropriate sizing level and call the corresponding function
    if level === :L2
        return alternative_sizing_L2(Q_hourly_to_three_pulses(Q), xy, rb, D, ks, Cs, s, ro, ri,
            kg, kp, kf, cf, ρf, μf, V, T0, Tlim; tp = tp, ny = ny, Tlim_at = Tlim_at)
    elseif level === :L3
        return alternative_sizing_L3(Q_hourly_to_monthly(Q), xy, rb, D, ks, Cs, s, ro, ri, kg, kp,
            kf, cf, ρf, μf, V, T0, Tlim; tp = tp, ny = ny, Tlim_at = Tlim_at)
    elseif level === :L4
        return alternative_sizing_L4(Q, xy, rb, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V,
            T0, Tlim; ny = ny, Tlim_at = Tlim_at)
    else
        throw(ArgumentError("level must be :L2, :L3 or :L4, got :$level"))
    end
end
