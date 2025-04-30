using .GHEModels

function alternative_sizing(Q₀::Vector{A}, ks::A, Cs::A, rb::A, D, Rb::A, xy::Matrix{A}, 
    T, method::String) where A<:Real
    """
    alternative_sizing(Q, ks, Cs, rb, D, Rb, xy, T, method)

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
    Outputs:
        - H: Borehole length for the lower and higher limit operating temperature (2x1) [m]
        - Out: Structure of various outputs parameters
    Reference:
        Ahmadfard, Mohammadamin, and Michel Bernier. 2019. “A Review of Vertical Ground Heat
        Exchanger Sizing Tools Including an Inter-Model Comparison.” Renewable and Sustainable 
        Energy Reviews 110:247–65. doi: 10.1016/j.rser.2019.04.045.
    """
    # Analyse the ground load for the method
    Q = Q_analysis(Q₀)

    if method == "L2"
        # Set heat loads for L2
        QL2 = [Q.y Q.y; Q.mₗ Q.mₕ; Q.hₗ Q.hₕ]
        # Perform sizing
        H = alternative_sizing_L2(QL2, ks, Cs, rb, D, Rb, xy, T)
    elseif method == "L3"
        # Set heat loads for L3
        QL3 = [Q.m̄ₐ, Q.m̄ₗ, Q.m̄ₕ]
        # Perform sizing
        H = alternative_sizing_L3(QL3, ks, Cs, rb, D, Rb, xy, T)
    elseif method == "L4"
        # Perform sizing
        H = alternative_sizing_L4(Q.hr, ks, Cs, rb, D, Rb, xy, T)
    elseif method == "all"
        H[1] = alternative_sizing_L2(Q, ks, Cs, rb, D, Rb, xy, T)
        H[2] = alternative_sizing_L3(Q, ks, Cs, rb, D, Rb, xy, T)
        H[3] = alternative_sizing_L4(Q, ks, Cs, rb, D, Rb, xy, T)
    else
        error("The method inputs can be 'L2', 'L3', 'L4' or 'all'")
    end
    return H
end

function alternative_sizing_L2(Q::Matrix{A}, ks::A, Cs::A, rb::A, D::A, Rb::A, xy::Matrix{A}, 
    T) where A<:Real
    """
        alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)

    Function that performs an 3 pulse (L2) sizing equation based on the alternative sizing equation.
    TODO: Add details for inputs and outputs
    """
    # 1. Define time array
    t = (y = 10 * 365 * 24 * 3600,
        m = round((((7 * 31) + (4 * 30) + 28) / 12) * 24 * 3600),
        h = 6 * 3600)

    #2. Iterate to find borehole length that covers the ground thermal loads
    Hᵢ = [0., 150.,]
    ii = 2
    Ltmp = Vector{Float64}(undef, 2)

    while abs(Hᵢ[ii] - Hᵢ[ii - 1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function

        g = bloc_matrix([t.h, t.m+t.h, t.y+t.m+t.h], ks, Cs, rb, Hᵢ[ii], D, xy)
        #g = successive_flux([t.h, t.m+t.h, t.y+t.m+t.h], k.s, Cs, r.b, Hi(ii), D, xy)*2*pi*k.s*nb;

        # Evaluate ground resistances R_gy, R_gm and R_gh
        Rgy = (g[3] - g[2]) / (2 * π * ks)
        Rgm = (g[2] - g[1]) / (2 * π * ks)
        Rgh = g[1] / (2 * pi * ks)

        # Evaluate sizing length
        Ltmp[1] = ((Q[1,1] * Rgy) + (Q[2,1] * Rgm) + (Q[3,1] * (Rgh + Rb))) / (T.L - T.g)
        Ltmp[2] = ((Q[1,2] * Rgy) + (Q[2,2] * Rgm) + (Q[3,2] * (Rgh + Rb))) / (T.H - T.g)
        Li = maximum(Ltmp)
        push!(Hᵢ, Li / size(xy, 2))
        print(Hᵢ)

        # Iteration output
        ii = ii + 1
        println("H = " * string(round(Hᵢ[ii]; digits=1)) * " m | tol = " *
                string(round(abs(Hᵢ[ii] - Hᵢ[ii - 1]); digits=1)) * " m")
    end
    return Hᵢ
end

function alternative_sizing_L3(Q, ks, Cs, rb, D, Rb, xy, T)
    """
    alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)

    Function that performs a monthly (L3) sizing equation based on the alternative sizing equation.
    """
    # 1. Define time array
    t = 3600:3600:(3600 * 24 * 365 * 10)
end

function alternative_sizing_L4(Q, k, C, r, D, Rb, xy, T)
    """
    alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)

    Function that performs an hourly (L4) sizing equation based on the alternative sizing equation.
"""
    # 1. Define time array
    t = 3600:3600:(3600 * 24 * 365 * 10)
end