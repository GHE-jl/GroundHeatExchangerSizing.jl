include("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

"""
    alternative_sizing_3pulses(Q, ks, Cs, rb, D, Rb, xy, T, method, t_peak, n_year)

Function that performs the alternative ASHRAE sizing equation to size a ground heat exchanger
(GHE) for either an three pulse (L2), monthly loads (L3) or hourly loads (L4)) approach. The direct
model should be based on the finite line source (FLS) model.
# Inputs:
    - `Q`: Ground heating power profile for one year (8760x1) [W]
    - `ks`: Soil thermal conductivity [W/mK]
    - `Cs`: Soil volumetric heat capacity [J/m³K]
    - `rb`: Borehole radius [m]
    - `D`: Borehole buried depth [m]
    - `Rb`: Borehole thermal resistance [mK/W]
    - `xy`: Borefield geometry (nbx2) [m]. nb is the number of boreholes. Each row has the x and y
        position of a borehole. If only one borehole is present, fill with [0,0].
    - `T₀`: Undisturbed ground temperature (1x1) [degC]
    - `Tlim`: Lower and higher operating fluid temperature limits (2x1) [degC]. E.g., [0.0, 30.0].
    - `tp`: Time of the peak load for L2 and L3 (default 6 hours) [hour]
    - `ny`: Number of years of the simulation (default 10 years) [year]
# Outputs:
    - `maximum(H)`: Longer borehole length obtained by sizing [m]
    - `out`: Structure of outputs parameters from the optimization
    - `H`: (Optional) Pair of length values for both temperature limits (2x1) [m]
# Reference:
    - Ahmadfard, M., & Bernier, M. (2018). Modifications to ASHRAE’s sizing method for vertical 
        ground heat exchangers. Science and Technology for the Built Environment, 24(7), 803–817. 
        https://doi.org/10.1080/23744731.2018.1423816
    - Ahmadfard, Mohammadamin, and Michel Bernier. 2019. “A Review of Vertical Ground Heat
        Exchanger Sizing Tools Including an Inter-Model Comparison.” Renewable and Sustainable 
        Energy Reviews 110:247–65. doi: 10.1016/j.rser.2019.04.045.
"""
function alternative_sizing_3pulses(Q::AbstractMatrix{<:Real}, ks::Real, Cs::Real, rb::Real, 
    D::Real, Rb::Real, xy::AbstractMatrix{<:Real}, T₀::Real, Tlim::AbstractVector{<:Real},
    tp::Real = 6.0, ny::Real = 10.0)
    # 0. Ensure type stability
    T = promote_type(eltype(Q), typeof(ks), typeof(Cs), typeof(ks), typeof(ks), typeof(ks), typeof(ks), typeof(ks), typeof(ks), typeof(ks), typeof(ks))

    # 1. Define time array
    t = (y = n_year * 365 * 24 * 3600.0,
        m = round((((7 * 31) + (4 * 30) + 28) / 12) * 24 * 3600.0),
        h = t_peak * 3600.0)

    #2. Iterate to find borehole length that covers the ground thermal loads
    Hᵢ = [0.0, 150.0]
    ii = 2
    Ltmp = Vector{Float64}(undef, 2)
    nb = size(xy, 1)

    while abs(Hᵢ[ii] - Hᵢ[ii - 1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function
        g = bloc_matrix([t.h, t.m + t.h, t.y + t.m + t.h], ks, Cs, rb, Hᵢ[ii], D, xy)
        #g = successive_flux([t.h, t.m+t.h, t.y+t.m+t.h], k.s, Cs, r.b, Hi(ii), D, xy)*2*pi*k.s*nb;

        # Evaluate ground resistances R_gy, R_gm and R_gh
        Rgy = (g[3] - g[2]) * nb
        Rgm = (g[2] - g[1]) * nb
        Rgh = g[1] * nb

        # Evaluate sizing length
        Ltmp[1] = ((Q[1, 1] * Rgy) + (Q[2, 1] * Rgm) + (Q[3, 1] * (Rgh + Rb))) / (T.L - T.g)
        Ltmp[2] = ((Q[1, 2] * Rgy) + (Q[2, 2] * Rgm) + (Q[3, 2] * (Rgh + Rb))) / (T.H - T.g)
        Li = maximum(Ltmp)
        push!(Hᵢ, Li / size(xy, 1))

        # Iteration output
        #println("H = " * string(round(Hᵢ[ii]; digits = 1)) * " m | tol = " *
        #        string(round(abs(Hᵢ[ii + 1] - Hᵢ[ii]); digits = 2)) * " m")
        ii += 1
    end
    return Hᵢ[end]
end

"""
    alternative_sizing(Q, ks, Cs, rb, D, Rb, xy, T, method, t_peak, n_year)

Function that performs the alternative ASHRAE sizing equation to size a ground heat exchanger
(GHE) for either an three pulse (L2), monthly loads (L3) or hourly loads (L4)) approach. The direct
model should be based on the finite line source (FLS) model.
# Inputs:
    - `Q`: Ground heating power profile for one year (8760x1) [W]
    - `ks`: Soil thermal conductivity [W/mK]
    - `Cs`: Soil volumetric heat capacity [J/m³K]
    - `rb`: Borehole radius [m]
    - `D`: Borehole buried depth [m]
    - `Rb`: Borehole thermal resistance [mK/W]
    - `xy`: Borefield geometry (nbx2) [m]. nb is the number of boreholes. Each row has the x and y
        position of a borehole. If only one borehole is present, fill with [0,0].
    - `T₀`: Undisturbed ground temperature (1x1) [degC]
    - `T.L`: Lower and higher operating fluid temperature limits (2x1) [degC]. E.g., [0.0, 30.0].
    - `tp`: Time of the peak load for L2 and L3 (default 6 hours) [hour]
    - `ny`: Number of years of the simulation (default 10 years) [year]
# Outputs:
    - `maximum(H)`: Longer borehole length obtained by sizing [m]
    - `out`: Structure of outputs parameters from the optimization
    - `H`: (Optional) Pair of length values for both temperature limits (2x1) [m]
# Reference:
    - Ahmadfard, M., & Bernier, M. (2018). Modifications to ASHRAE’s sizing method for vertical 
        ground heat exchangers. Science and Technology for the Built Environment, 24(7), 803–817. 
        https://doi.org/10.1080/23744731.2018.1423816
    - Ahmadfard, Mohammadamin, and Michel Bernier. 2019. “A Review of Vertical Ground Heat
        Exchanger Sizing Tools Including an Inter-Model Comparison.” Renewable and Sustainable 
        Energy Reviews 110:247–65. doi: 10.1016/j.rser.2019.04.045.
"""
function alternative_sizing_L3(Q::Matrix{F}, ks::F, Cs::F, rb::F, D::F, Rb::F, xy::Matrix{F},
    T, t_peak = 6.0, n_year = 10.0) where {F <: AbstractFloat}
    # 0. Check optional input
    isa(t_peak, Float64) ? t_peak : Float64(n_year)
    isa(n_year, Float64) ? n_year : Float64(n_year)

    # 1. Define time array
    t = collect(3600.0:3600.0:(3600.0 * 24 * 365 * n_year))

    # 2. Set the loads for the L3 method
    Q = Q_L3_resample(Q, t_peak, n_year)
    f = diff([[0.0 0.0]; Q], dims=1)

    # 3. Iterate to find borehole length that covers the ground thermal loads
    Hᵢ = [0.0, 150.0]
    ii = 2
    Ltmp = Vector{Float64}(undef, 2)
    nb = size(xy, 1)

    while abs(Hᵢ[ii] - Hᵢ[ii - 1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function
        g = bloc_matrix(t, ks, Cs, rb, Hᵢ[ii], D, xy)

        # Evaluate sizing length
        Ltmp[1] = minimum(convolution(f[:, 1], g * nb) + (Q[:, 1] * Rb)) / (T.L - T.g)
        Ltmp[2] = maximum(convolution(f[:, 2], g * nb) + (Q[:, 2] * Rb)) / (T.H - T.g)
        Li = maximum(Ltmp)
        push!(Hᵢ, Li / size(xy, 1))

        # Iteration output
        #println("H = " * string(round(Hᵢ[ii]; digits = 1)) * " m | tol = " *
        #        string(round(abs(Hᵢ[ii + 1] - Hᵢ[ii]); digits = 2)) * " m")
        ii += 1
    end
    return Hᵢ[end]
end

"""
    alternative_sizing(Q, ks, Cs, rb, D, Rb, xy, T, method, t_peak, n_year)

Function that performs the alternative ASHRAE sizing equation to size a ground heat exchanger
(GHE) for either an three pulse (L2), monthly loads (L3) or hourly loads (L4)) approach. The direct
model should be based on the finite line source (FLS) model.
# Inputs:
    - `Q`: Ground heating power profile for one year (8760x1) [W]
    - `ks`: Soil thermal conductivity [W/mK]
    - `Cs`: Soil volumetric heat capacity [J/m³K]
    - `rb`: Borehole radius [m]
    - `D`: Borehole buried depth [m]
    - `Rb`: Borehole thermal resistance [mK/W]
    - `xy`: Borefield geometry (nbx2) [m]. nb is the number of boreholes. Each row has the x and y
        position of a borehole. If only one borehole is present, fill with [0,0].
    - `T₀`: Undisturbed ground temperature (1x1) [degC]
    - `T.L`: Lower and higher operating fluid temperature limits (2x1) [degC]. E.g., [0.0, 30.0].
    - `tp`: Time of the peak load for L2 and L3 (default 6 hours) [hour]
    - `ny`: Number of years of the simulation (default 10 years) [year]
# Outputs:
    - `maximum(H)`: Longer borehole length obtained by sizing [m]
    - `out`: Structure of outputs parameters from the optimization
    - `H`: (Optional) Pair of length values for both temperature limits (2x1) [m]
# Reference:
    - Ahmadfard, M., & Bernier, M. (2018). Modifications to ASHRAE’s sizing method for vertical 
        ground heat exchangers. Science and Technology for the Built Environment, 24(7), 803–817. 
        https://doi.org/10.1080/23744731.2018.1423816
    - Ahmadfard, Mohammadamin, and Michel Bernier. 2019. “A Review of Vertical Ground Heat
        Exchanger Sizing Tools Including an Inter-Model Comparison.” Renewable and Sustainable 
        Energy Reviews 110:247–65. doi: 10.1016/j.rser.2019.04.045.
"""
function alternative_sizing_L4(Q₀::AbstractVector{<:Real}, V::Real, ks::Real, kg::Real, kp::Real,
    kf::Real, Cs::Real, rb::Real, ro::Real, ri::Real, s::Real, D::Real, cf::Real, ρf::Real, 
    μf::Real, xy::AbstractArray{<:Real}, T, n_year::Real = 10.0)
    # 0. Check optional input
    isa(n_year, Float64) ? n_year : Float64(n_year)

    # 1. Define time array
    t = 3600.0:3600.0:(3600.0 * 24 * 365 * n_year)

    # 2. Precompute some values for the thermal resistance
    Rf = R_f(V / (size(xy, 1) * π * ri^2), kf, ri, cf, ρf, μf)
    Rp = R_p(kp, ro, ri)
    Rb = R_b_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)
    Ra = R_a_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    # 3. Iterate to find borehole length that covers the ground thermal loads
    L_tmp = Vector{Float64}(undef, 2)               # Preallocation
    g = Vector{Float64}(undef, 8760*Int(n_year))    # Preallocation
    tmp = similar(g)                                # Preallocation
    Rbₑ = zero(Float64)                             # Preallocation
    Q = repeat(Q₀, Int(n_year))                     # Repeated thermal loads
    f = diff([0.0; Q])                              # Impulse function
    Hᵢ = [0.0, 150.0]                               # Initialize the borehole length
    nb = size(xy, 1)                                # Number of boreholes
    ii = 2                                          # Start the iterative counter

    while abs(Hᵢ[ii] - Hᵢ[ii - 1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function for impulse in W/m (so in units °Cm/W)
        #g = bloc_matrix(t, ks, Cs, rb, Hᵢ[ii], D, xy)
        g = g_model(t, ks, Cs, rb, Hᵢ[ii], D, xy, model = "fls")

        # Evaluate the effective borehole thermal resistance
        Rbₑ = R_bₑ(V, cf, ρf, Hᵢ[ii], Rb, Ra)

        # Evaluate sizing length
        #tmp = convolution(f, g * 2 * π * ks * size(xy, 1)) .+ (Q * Rbₑ) # Dimensionless g
        tmp = convolution(f, g * nb) .+ (Q * Rbₑ)    # g in °Cm/W
        L_tmp[1] = minimum(tmp) / (T.L - T.g)
        L_tmp[2] = maximum(tmp) / (T.H - T.g)
        Li = maximum(L_tmp)
        push!(Hᵢ, Li / size(xy, 1))

        # Iteration output
        #=println("H = " * string(round(Hᵢ[ii]; digits = 1)) * " m | tol = " *
                string(round(abs(Hᵢ[ii + 1] - Hᵢ[ii]); digits = 2)) * " m")=#
        ii += 1
    end
    return Hᵢ[end]
end