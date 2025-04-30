"""
Script performing sizing using the alternative ASHRAE sizing equation.
"""

using BenchmarkTools
using Plots

# Add GHE models
includet("../../GHEModels.jl/src/GHEModels.jl")

# Add the sizing equations
includet("../src/GHESizings.jl")
using .GHESizings

# Add the cases from Ahmadfard and Bernier (2019)
includet("Ahmadfard_Cases.jl")

# 1. Select case
Q₀, n, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_1()
#Q₀, n, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_2()
#Q₀, n, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_3()
#Q₀, n, nb, B, D, r, s, ρ, C, μ, k, T, R, V = case_4()

# 2. Analyse the heat load profile
Q = Q_analysis(Q₀)

# 3. Evaluate the borehole thermal resistance
Rb = Rb_first_order_multipole(
    V / (nb * π * r.i^2), k.s, k.g, k.p, k.f, C.f / ρ.f, ρ.f, μ, r.b, r.i, r.o, s)

# 2. Call the alternative sizing equation
QL2 = [Q.y Q.y; Q.mₗ Q.mₕ; Q.hₗ Q.hₕ]
xy = B * hcat([[i, j] for i in 1:n.x for j in 1:n.y]...)'.-B
H = alternative_sizing_L2(QL2, k.s, C.s, r.b, D, Rb, xy, T)