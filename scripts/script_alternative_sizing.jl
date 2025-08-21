"""
Script performing sizing using the alternative ASHRAE sizing equation.
"""

using BenchmarkTools
#using CairoMakie

includet("../../GHEModels.jl/src/GHEModels.jl")
includet("../src/GHESizings.jl")
includet("Ahmadfard_Cases.jl")
using .GHEModels, .GHESizings

# 1. Select case
#Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_1()
#Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_2()
#Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_3()
Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_4()

Q₀ = Float64.(Q₀)

# 2. Analyse the heat load profile
Q = Q_analysis(Q₀)
t_peak = 6.0
n_year = 10.0

# 3. Evaluate the borehole thermal resistance
Rf = R_f(V / (size(xy, 1) * π * r.i^2), k.f, r.i, C.f / ρ.f, ρ.f, μ)
Rp = R_p(k.p, r.o, r.i)
Rb = R_b_first_order_multipole(k.s, k.g, r.b, r.o, s, Rp, Rf)
Ra = R_a_first_order_multipole(k.s, k.g, r.b, r.o, s, Rp, Rf)
Rbₑ = R_bₑ(V, C.f / ρ.f, ρ.f, 100.0, Rb, Ra)

# 4. Call the L2 alternative sizing equation
QL2 = [Q.y Q.y; Q.mₗ Q.mₕ; Q.hₗ Q.hₕ]
@time H2 = alternative_sizing_L2(QL2, k.s, C.s, r.b, D, Rb, xy, T)

# 5. Call the L3 alternative sizing equation
QL3 = [Q.m̄ₐ Q.m̄ₗ Q.m̄ₕ]
@time H3 = alternative_sizing_L3(QL3, k.s, C.s, r.b, D, Rb, xy, T)

# 6. Call the L4 alternative sizing equation
@time H4 = alternative_sizing_L4(Q₀, V, k.s, k.g, k.p, k.f, C.s, r.b, r.o, r.i, s, D, C.f / ρ.f, 
    ρ.f, μ, xy, T, n_year)

# 7. Call all method with the master function
#H = alternative_sizing(Q.hr, k.s, C.s, r.b, D, Rb, xy, T)