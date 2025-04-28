using GHEModels

function alternative_sizing(Q, k, C, r, D, Rb, xy, T, method)
    """
    alternative_sizing(Q, k, C, r, D, Rb, xy, T, method)
    
    Function that performs the alternative ASHRAE sizing equation to size a ground heat exchanger
    (GHE) for either an L2 (three pulse), L3 (monthly loads) or L4 (hourly loads) approach
    (Ahmadfard and Bernier, 2019). The direct model should be based on the finite line source (FLS)
    model.
    Inputs:
        - Q: Ground heating power profile for one year (W) [8760x1]
        - k: Structure of thermal conductivity (W/mK):
            - k.s: Soil [1x1]
            - k.g: Grout [1x1]
            - k.p: Pipe [1x1]
        - Cs: Soil volumetric specific heat (J/m^3K) [1x1]
        - r: Structure of radius (m):
            - r.b: Borehole [1x1]
            - r.i: Pipe inlet [1x1]
            - r.o: Pipe outlet [1x1]
        - D: Borehole buried depth (m) [1x1]
        - Rb: Borehole thermal resistance (mK/W) [1x1]
        - xy: Borefield geometry (m) [nbx2]. nb is the number of boreholes. Each row has the x and y
            position of a borehole. If only one borehole is present, fill with [0,0].
        - Tg: Undisturbed ground temperature (degC) [1x1]
        - Tlim: Lower and higher limit operating average fluid temperature (degC) [2x1]. The first
            input is the lower limit and the second is the higher. E.g., [0, 35]
        - Method: Either 'L2', 'L3' or 'L4' for the 3 available sizing equation.
    Outputs:
        - H: Borehole length for the lower and higher limit operating temperature (m) [2x1]
        - Out: Structure of various outputs parameters
    Reference:
        Ahmadfard, Mohammadamin, and Michel Bernier. 2019. “A Review of Vertical Ground Heat
        Exchanger Sizing Tools Including an Inter-Model Comparison.” Renewable and Sustainable 
        Energy Reviews 110:247–65. doi: 10.1016/j.rser.2019.04.045.
    """

    if method == "L2"
        H = alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)
    elseif method == "L3"
        H = alternative_sizing_L3(Q, k, C, r, D, Rb, xy, T)
    elseif method == "L4"
        H = alternative_sizing_L4(Q, k, C, r, D, Rb, xy, T)
    elseif method == "all"
        H[1] = alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)
        H[2] = alternative_sizing_L3(Q, k, C, r, D, Rb, xy, T)
        H[3] = alternative_sizing_L4(Q, k, C, r, D, Rb, xy, T)
    else
        error("The method inputs can be 'L2', 'L3', 'L4' or 'all'")
    end
    return H
end

function alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)
    """
        alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)
    
    Function that performs an 3 pulse (L2) sizing equation based on the alternative sizing equation.
    """
    # 1. Define time array
    t = (y = 10*365*24*3600,
        m = round((((7*31)+(4*30)+28)/12)*24*3600),
        h = 6*3600)
    
    #2. Iterate to find borehole length that covers the ground thermal loads
    Hᵢ = [0, 150]
    ii = 2

    while abs(Hi[ii] - Hi[ii-1]) > 0.01 && ii < 20  # 0.01 m and 20 iterations for convergence
        # Evaluate the g-function
        g = bloc_matrix(params::GHEParam, xy::Matrix{T})
        g = successive_flux([t.h, t.m+t.h, t.y+t.m+t.h],
            k.s, Cs, r.b, Hi(ii), D, xy)*2*pi*k.s*nb;
    
        # Evaluate ground resistances R_gy, R_gm and R_gh
        R.gy = (g(3)-g(2))/(2*pi*k.s);
        R.gm = (g(2)-g(1))/(2*pi*k.s);
        R.gh = g(1)/(2*pi*k.s);
    
        # Evaluate sizing length
        Ltmp(1) = (Q.y*R.gy + Q.m_c*R.gm + Q.h_c*(R.gh+Rb))/(Tlim(1)-Tg);
        Ltmp(2) = (Q.y*R.gy + Q.m_h*R.gm + Q.h_h*(R.gh+Rb))/(Tlim(2)-Tg);
        Li = max(Ltmp);
        Hi(ii+1) = Li/nb;
    
        # Iteration output
        ii = ii+1;
        println("H = "*round(Hᵢ[ii], 1)*" m | tol = "*round(abs(Hi[ii] - Hi[ii-1]), 1)*" m")
    end

end


function alternative_sizing_L3(Q, k, C, r, D, Rb, xy, T)
    """
    alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)

    Function that performs a monthly (L3) sizing equation based on the alternative sizing equation.
    """
    # 1. Define time array
    t = 3600:3600:3600*24*365*10;

end

function alternative_sizing_L4(Q, k, C, r, D, Rb, xy, T)
    """
    alternative_sizing_L2(Q, k, C, r, D, Rb, xy, T)

    Function that performs an hourly (L4) sizing equation based on the alternative sizing equation.
"""
    # 1. Define time array
    t = 3600:3600:3600*24*365*10;

end