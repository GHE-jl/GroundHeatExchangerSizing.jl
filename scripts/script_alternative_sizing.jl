"""
Script performing sizing using the alternative ASHRAE sizing equation.
"""

using BenchmarkTools
using Plots

# Add GHE models
includet("../../GHEModels.jl/src/GHEModels.jl")
using .GHEModels

# Add the sizing equations
includet("../src/GHESizings.jl")
using .GHESizings

# Add the cases from Ahmadfard and Bernier (2019)
includet("Ahmadfard_Cases.jl")

# 1. Select case
Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_1()
#Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_2()
#Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_3()
#Q₀, n, xy, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_4()

# 2. Analyse the heat load profile
Q = Q_analysis(Q₀)
t_peak = 6.0
n_year = 10.0

# 3. Evaluate the borehole thermal resistance
Rb = Rb_first_order_multipole(
    V / (nb * π * r.i^2), k.s, k.g, k.p, k.f, C.f / ρ.f, ρ.f, μ, r.b, r.i, r.o, s)

# 4. Call the L2 alternative sizing equation
#QL2 = [Q.y Q.y; Q.mₗ Q.mₕ; Q.hₗ Q.hₕ]
#@time H2 = alternative_sizing_L2(QL2, k.s, C.s, r.b, D, Rb, xy, T)

# 5. Call the L3 alternative sizing equation
#QL3 = [Q.m̄ₐ Q.m̄ₗ Q.m̄ₕ]
#@time H3 = alternative_sizing_L3(QL3, k.s, C.s, r.b, D, Rb, xy, T)

# 6. Call the L4 alternative sizing equation
#@time H4 = alternative_sizing_L4(Q.hr, k.s, C.s, r.b, D, Rb, xy, T)

# 7. Call all method with the master function
H = alternative_sizing(Q.hr, k.s, C.s, r.b, D, Rb, xy, T)