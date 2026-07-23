# Borehole-outlet transfer-function sizing (Dion & Pasquier 2025) on the Ahmadfard & Bernier (2019)
# test cases (as adapted in Dion & Pasquier 2025, Table 1). Sizes each case with the L2, L3 and L4
# levels and compares the lengths with the alternative ASHRAE equation.

using Pkg; Pkg.activate(@__DIR__)
Pkg.instantiate()

using GroundHeatExchangerSizing
using CairoMakie
using BenchmarkTools
using Logging

# The short-term ANN warns on every evaluation for the test-case pipe/fluid properties that fall
# outside its training ranges (see the note in the docs). Those inputs do not vary during the
# length search, so silence Warn-level logging for this run instead of printing thousands of copies.
disable_logging(Logging.Warn)

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
        H_out[ci, li] = outlet_sizing(Q, p.xy, p.r.b, p.D, p.k.s, p.C.s, p.s, p.r.o, p.r.i, p.k.g,
            p.C.g, p.k.p, p.C.p, p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ, V, p.T.g, [p.T.L, p.T.H];
            level = lvl)
        H_alt[ci, li] = alternative_sizing(Q, p.xy, p.r.b, p.D, p.k.s, p.C.s, p.s, p.r.o, p.r.i,
            p.k.g, p.k.p, p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ, V, p.T.g, [p.T.L, p.T.H]; level = lvl)
        println("Case $c  $lvl : outlet H = $(round(H_out[ci, li]; digits = 1)) m | " *
                "alternative H = $(round(H_alt[ci, li]; digits = 1)) m")
    end
end

# Comparison against Dion & Pasquier (2025), Table 2 — governing (larger of low/high) lengths [m].
H_ref_gfun = [59.2 58.0 56.6; 86.3 88.7 91.8; 94.3 117.9 118.8; 111.4 111.3 111.1]
H_ref_outlet = [53.0 53.9 50.2; 83.0 85.0 88.9; 90.0 113.6 114.8; 106.0 106.5 104.2]

relative_error(H, H_ref) = (H .- H_ref) ./ H_ref .* 100

RE_alt = relative_error(H_alt, H_ref_gfun)
RE_out = relative_error(H_out, H_ref_outlet)

println("\nRelative error vs Dion & Pasquier (2025), Table 2 [%]")
println("Case  Level  alternative vs g-functions   outlet vs transfer functions")
for (ci, c) in enumerate(cases), (li, lvl) in enumerate(levels)
    println("  $c    $lvl        $(round(RE_alt[ci, li]; digits = 1))" *
            "                          $(round(RE_out[ci, li]; digits = 1))")
end
println("Mean |RE| — alternative: $(round(sum(abs, RE_alt) / length(RE_alt); digits = 2))% | " *
        "outlet: $(round(sum(abs, RE_out) / length(RE_out); digits = 2))%")

# Outlet sizing time per level (BenchmarkTools), on case 4 (5×5 field).
#=
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
=#

# Compare the two sizing families per level.
fig = Figure(size = (900, 320))
for (li, lvl) in enumerate(levels)
    ax = Axis(fig[1, li], xlabel = "Test case", ylabel = li == 1 ? "H [m]" : "",
        title = string(lvl), xticks = (collect(cases), string.(cases)))
    barplot!(ax, collect(cases) .- 0.2, H_alt[:, li]; width = 0.4, label = "alternative")
    barplot!(ax, collect(cases) .+ 0.2, H_out[:, li]; width = 0.4, label = "outlet")
    li == 3 && axislegend(ax; position = :lt)
end
display(fig)