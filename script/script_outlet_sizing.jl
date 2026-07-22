# Borehole-outlet transfer-function sizing (Dion & Pasquier 2025) on the Ahmadfard & Bernier (2019)
# test cases (as adapted in Dion & Pasquier 2025, Table 1). Sizes each case with the L2, L3 and L4
# levels and compares the lengths with the alternative ASHRAE equation.
#
# Run from the package root with the script environment:
#   julia --project=script -e 'using Pkg; Pkg.instantiate()'
#   julia --project=script script/script_outlet_sizing.jl

using GroundHeatExchangerSizing
using CairoMakie
using BenchmarkTools

include("Ahmadfard_cases.jl")

levels = (:L2, :L3, :L4)
cases = 1:4
H_out = zeros(length(cases), length(levels))
H_alt = zeros(length(cases), length(levels))

for (ci, c) in enumerate(cases)
    p = ahmadfard_cases(c)
    Q = Float64.(p.Q)
    V = p.V / p.nb
    for (li, lvl) in enumerate(levels)
        r_out = outlet_sizing(Q, p.xy, p.r.b, p.D, p.k.s, p.C.s, p.s, p.r.o, p.r.i, p.k.g, p.C.g,
            p.k.p, p.C.p, p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ, V, p.T.g, [p.T.L, p.T.H]; level = lvl)
        r_alt = alternative_sizing(Q, p.xy, p.r.b, p.D, p.k.s, p.C.s, p.s, p.r.o, p.r.i, p.k.g,
            p.k.p, p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ, V, p.T.g, [p.T.L, p.T.H]; level = lvl)
        H_out[ci, li] = r_out.H
        H_alt[ci, li] = r_alt.H
        println("Case $c  $lvl : outlet H = $(round(r_out.H; digits = 1)) m | " *
                "alternative H = $(round(r_alt.H; digits = 1)) m")
    end
end

# Outlet sizing time per level (BenchmarkTools), on case 4 (5×5 field). The outlet method solves an
# optimisation per limit, so it is heavier than the alternative fixed-point; L3/L4 stay tractable
# because the transfer function is sub-sampled and PCHIP-interpolated (interp = true).
let p = ahmadfard_cases(4)
    Q = Float64.(p.Q); xy = p.xy; V = p.V / p.nb
    rb, D, ks, Cs, s = p.r.b, p.D, p.k.s, p.C.s, p.s
    ro, ri, kg, Cg, kp, Cp = p.r.o, p.r.i, p.k.g, p.C.g, p.k.p, p.C.p
    kf, cf, ρf, μf = p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ
    T0, Tlim = p.T.g, [p.T.L, p.T.H]
    println("\nOutlet sizing time per level (case 4):")
    for lvl in levels
        print("  ", lvl, ": ")
        @btime outlet_sizing($Q, $xy, $rb, $D, $ks, $Cs, $s, $ro, $ri, $kg, $Cg, $kp, $Cp,
            $kf, $cf, $ρf, $μf, $V, $T0, $Tlim; level = $lvl)
    end
end

# Compare the two sizing families per level.
fig = Figure(size = (900, 320))
for (li, lvl) in enumerate(levels)
    ax = Axis(fig[1, li], xlabel = "Test case", ylabel = li == 1 ? "H [m]" : "",
        title = string(lvl), xticks = (collect(cases), string.(cases)))
    barplot!(ax, collect(cases) .- 0.2, H_alt[:, li]; width = 0.4, label = "alternative")
    barplot!(ax, collect(cases) .+ 0.2, H_out[:, li]; width = 0.4, label = "outlet")
    li == 3 && axislegend(ax; position = :lt)
end

mkpath(joinpath(@__DIR__, "figures"))
save(joinpath(@__DIR__, "figures", "outlet_sizing.png"), fig)
println("Saved figures/outlet_sizing.png")
