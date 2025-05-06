#using .GHEModels

function alternative_sizing(Q₀::Vector{A}, ks::A, Cs::A, rb::A, D, Rb::A, xy::Matrix{A},
        T, method = "all", t_peak = nothing, n_year = nothing) where {A <: Real}
    """
    alternative_sizing(Q, ks, Cs, rb, D, Rb, xy, T, method, t_peak, n_year)

    Function that performs the alternative ASHRAE sizing equation to size a ground heat exchanger
    (GHE) for either an L2 (three pulse), L3 (monthly loads) or L4 (hourly loads) approach
    (Ahmadfard and Bernier, 2019). The direct model should be based on the finite line source (FLS)
    model.
    Inputs:
        - Q: Ground heating power profile for one year (8760x1) [W]
        - ks: Soil thermal conductivity (1x1) [W/mK]
        - Cs: Soil volumetric specific heat (1x1) [J/m³K]
        - rb: Borehole radius (1x1) [m]
        - D: Borehole buried depth (1x1) [m]
        - Rb: Borehole thermal resistance (1x1) [mK/W]
        - xy: Borefield geometry (nbx2) [m]. nb is the number of boreholes. Each row has the x and y
            position of a borehole. If only one borehole is present, fill with [0,0].
        - T: Temperature [degC]
            - T.g: Undisturbed ground temperature (1x1) [degC]
            - T.L: Lower operating fluid temperature limit (1x1) [degC]. E.g., 0.0
            - T.H: Higher operating fluid temperature limit (1x1) [degC]. E.g., 30.0
        - Method: Either 'L2', 'L3' or 'L4' for the 3 available sizing equation.
        - t_peak: (optionnal) Time of the peak load for L2 and L3 (default 6 hours) [hour]
        - n_year: (optionnal) Number of years of the simulation (default: 10 years) [year]
    Outputs:
        - H: Borehole length for the lower and higher limit operating temperature (2x1) [m]
        - Out: Structure of various outputs parameters
    Reference:
        Ahmadfard, Mohammadamin, and Michel Bernier. 2019. “A Review of Vertical Ground Heat
        Exchanger Sizing Tools Including an Inter-Model Comparison.” Renewable and Sustainable 
        Energy Reviews 110:247–65. doi: 10.1016/j.rser.2019.04.045.
    """
    # Set optional Inputs
    isnothing(t_peak) ? t_peak = 6.0 : t_peak
    isnothing(n_year) ? n_year = 10.0 : n_year

    # Analyse the ground load for the method
    Q = Q_analysis(Q₀)

    if method == "L2"
        QL2 = [Q.y Q.y; Q.mₗ Q.mₕ; Q.hₗ Q.hₕ]   # Set loads
        H = alternative_sizing_L2(QL2, ks, Cs, rb, D, Rb, xy, T, t_peak, n_year)
    elseif method == "L3"
        QL3 = [Q.m̄ₐ Q.m̄ₗ Q.m̄ₕ]                  # Set loads
        H = alternative_sizing_L3(QL3, ks, Cs, rb, D, Rb, xy, T, t_peak, n_year)
    elseif method == "L4"
        H = alternative_sizing_L4(Q.hr, ks, Cs, rb, D, Rb, xy, T, n_year)
    elseif method == "all"
        H = Vector{Float64}(undef, 3)           # Preallocating output
        QL2 = [Q.y Q.y; Q.mₗ Q.mₕ; Q.hₗ Q.hₕ]      # Loads for L2
        QL3 = [Q.m̄ₐ Q.m̄ₗ Q.m̄ₕ]                  # Loads for L3
        H[1] = alternative_sizing_L2(QL2, ks, Cs, rb, D, Rb, xy, T, t_peak, n_year)
        H[2] = alternative_sizing_L3(QL3, ks, Cs, rb, D, Rb, xy, T, t_peak, n_year)
        H[3] = alternative_sizing_L4(Q.hr, ks, Cs, rb, D, Rb, xy, T, n_year)
    else
        error("The method inputs can be 'L2', 'L3', 'L4' or 'all'")
    end
    return H
end

function alternative_sizing_L2(Q::Matrix{A}, ks::A, Cs::A, rb::A, D::A, Rb::A, xy::Matrix{A},
        T, t_peak = nothing, n_year = nothing) where {A <: Real}
    """
        alternative_sizing_L2(Q, ks, Cs, rb, D, Rb, xy, T, t_peak=6.0, n_year=10)

    Function that performs an 3 pulse (L2) sizing equation based on the alternative sizing equation.
    The heat load variable Q is a (3x2) matrix, with the first column for ground cooling loads, and
    the second column for ground heating load. The first row is yearly mean, the second row is the
    monthly average load when the peak load occurs, and the third row is the maximum load.
    """
    # 0. Set optional Inputs
    isnothing(t_peak) ? t_peak = 6.0 : t_peak
    isnothing(n_year) ? n_year = 10.0 : n_year

    # 1. Define time array
    t = (y = n_year * 365 * 24 * 3600.0,
        m = round((((7 * 31) + (4 * 30) + 28) / 12) * 24 * 3600.0),
        h = t_peak * 3600.0)

    #2. Iterate to find borehole length that covers the ground thermal loads
    Hᵢ = [0.0, 150.0]
    ii = 2
    Ltmp = Vector{Float64}(undef, 2)

    while abs(Hᵢ[ii] - Hᵢ[ii - 1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function
        g = bloc_matrix([t.h, t.m + t.h, t.y + t.m + t.h], ks, Cs, rb, Hᵢ[ii], D, xy)
        #g = successive_flux([t.h, t.m+t.h, t.y+t.m+t.h], k.s, Cs, r.b, Hi(ii), D, xy)*2*pi*k.s*nb;

        # Evaluate ground resistances R_gy, R_gm and R_gh
        Rgy = (g[3] - g[2]) / (2 * π * ks)
        Rgm = (g[2] - g[1]) / (2 * π * ks)
        Rgh = g[1] / (2 * pi * ks)

        # Evaluate sizing length
        Ltmp[1] = ((Q[1, 1] * Rgy) + (Q[2, 1] * Rgm) + (Q[3, 1] * (Rgh + Rb))) / (T.L - T.g)
        Ltmp[2] = ((Q[1, 2] * Rgy) + (Q[2, 2] * Rgm) + (Q[3, 2] * (Rgh + Rb))) / (T.H - T.g)
        Li = maximum(Ltmp)
        push!(Hᵢ, Li / size(xy, 1))

        # Iteration output
        println("H = " * string(round(Hᵢ[ii]; digits = 1)) * " m | tol = " *
                string(round(abs(Hᵢ[ii + 1] - Hᵢ[ii]); digits = 2)) * " m")
        ii += 1
    end
    return Hᵢ[end]
end

function alternative_sizing_L3(Q::Matrix{A}, ks::A, Cs::A, rb::A, D::A, Rb::A, xy::Matrix{A},
    T, t_peak = nothing, n_year = nothing) where {A <: Real}
    """
    alternative_sizing_L3(Q, ks, Cs, rb, D, Rb, xy, T; t_peak=6.0, n_year=10)

    Function that performs a monthly (L3) sizing equation based on the alternative sizing equation.
    The hear load variable Q is a matrix (12x3) formed of [Q.m̄ₐ, Q.m̄ₗ, Q.m̄ₕ] from the QLoads 
    structure, meaning the average monthly loads and the minimum and maximum laods.
    """
    # 0. Set optional Inputs
    isnothing(t_peak) ? t_peak = 6.0 : t_peak
    isnothing(n_year) ? n_year = 10.0 : n_year

    # 1. Define time array
    t = collect(3600.0:3600.0:(3600.0 * 24 * 365 * n_year))

    # 2. Set the loads for the L3 method
    Q = Q_L3_resample(Q, t_peak, n_year)
    f = diff([[0.0 0.0]; Q], dims=1)

    # 3. Iterate to find borehole length that covers the ground thermal loads
    Hᵢ = [0.0, 150.0]
    ii = 2
    Ltmp = Vector{Float64}(undef, 2)

    while abs(Hᵢ[ii] - Hᵢ[ii - 1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function
        g = bloc_matrix(t, ks, Cs, rb, Hᵢ[ii], D, xy)

        # Evaluate sizing length
        Ltmp[1] = minimum(convolution(f[:, 1], g/(2*π*ks)) + (Q[:, 1] * Rb)) / (T.L - T.g)
        Ltmp[2] = maximum(convolution(f[:, 2], g/(2*π*ks)) + (Q[:, 2] * Rb)) / (T.H - T.g)
        Li = maximum(Ltmp)
        push!(Hᵢ, Li / size(xy, 1))

        # Iteration output
        println("H = " * string(round(Hᵢ[ii]; digits = 1)) * " m | tol = " *
                string(round(abs(Hᵢ[ii + 1] - Hᵢ[ii]); digits = 2)) * " m")
        ii += 1
    end
    return Hᵢ[end]
end

function alternative_sizing_L4(Q₀::Vector{A}, ks::A, Cs::A, rb::A, D::A, Rb::A, xy::Matrix{A},
    T, n_year = nothing) where {A <: Real}
    """
    alternative_sizing_L4(Q₀, ks, Cs, rb, D, Rb, xy, T; n_year=10)

    Function that performs a hourly (L4) sizing equation based on the alternative sizing equation.
    """
    # 0. Set optional Inputs
    isnothing(n_year) ? n_year = 10.0 : n_year

    # 1. Define time array
    t = collect(3600.0:3600.0:(3600.0 * 24 * 365 * n_year))

    # 2. Iterate to find borehole length that covers the ground thermal loads
    L_tmp = Vector{Float64}(undef, 2)               # Preallocation
    g = Vector{Float64}(undef, 8760*Int(n_year))    # Preallocation
    tmp = similar(g)                                # Preallocation
    Q = repeat(Q₀, Int(n_year))
    f = diff([0.0; Q])
    Hᵢ = [0.0, 150.0]
    ii = 2

    while abs(Hᵢ[ii] - Hᵢ[ii - 1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function
        g = bloc_matrix(t, ks, Cs, rb, Hᵢ[ii], D, xy)

        # Evaluate sizing length
        tmp = convolution(f, g/(2*π*ks)) .+ (Q * Rb)
        L_tmp[1] = minimum(tmp) / (T.L - T.g)
        L_tmp[2] = maximum(tmp) / (T.H - T.g)
        Li = maximum(L_tmp)
        push!(Hᵢ, Li / size(xy, 1))

        # Iteration output
        println("H = " * string(round(Hᵢ[ii]; digits = 1)) * " m | tol = " *
                string(round(abs(Hᵢ[ii + 1] - Hᵢ[ii]); digits = 2)) * " m")
        ii += 1
    end
    return Hᵢ[end]
end