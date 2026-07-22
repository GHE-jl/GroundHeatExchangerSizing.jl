# Verification of the alternative ASHRAE sizing equation (g-function based) on the four Ahmadfard &
# Bernier (2019) test cases. Sizes each case with the L2, L3 and L4 fixed-point iterations and
# compares the borehole lengths against the g-function reference lengths of Dion & Pasquier (2025),
# Table 2 (the "g-functions" columns).

using Pkg; Pkg.activate(@__DIR__)
Pkg.instantiate()

using GroundHeatExchangerSizing
using CairoMakie
using BenchmarkTools

# Reference cases (hourly ground loads + all parameters) from Ahmadfard & Bernier (2019).
include("Ahmadfard_cases.jl")

levels = (:L2, :L3, :L4)
cases = 1:4

# Reference lengths [m] — Dion & Pasquier (2025), Table 2, g-function columns (L2, L3, L4).
H_ref = [59.2 58.0 56.6; 86.3 88.7 91.8; 94.3 117.9 118.8; 111.4 111.3 111.1]

H = zeros(length(cases), length(levels))

println("Case  level   H (m)   ref (m)   Δ (%)")
for (ci, c) in enumerate(cases)
    p = ahmadfard_cases(c)
    Q = Float64.(p.Q)
    V = p.V / p.nb                     # per-borehole loop flow [m³/s] (case stores total flow)
    for (li, lvl) in enumerate(levels)
        res = alternative_sizing(Q, p.xy, p.r.b, p.D, p.k.s, p.C.s, p.s, p.r.o, p.r.i, p.k.g,
            p.k.p, p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ, V, p.T.g, [p.T.L, p.T.H]; level = lvl)
        H[ci, li] = res.H
        Δ = 100 * (res.H - H_ref[ci, li]) / H_ref[ci, li]
        println(lpad(c, 3), lpad(string(lvl), 6), lpad(round(res.H; digits = 1), 8),
            lpad(round(H_ref[ci, li]; digits = 1), 9), lpad(round(Δ; digits = 1), 8))
    end
end

# Sizing time per level (BenchmarkTools), on case 4 (5×5 field). Interpolating locals ($...) keeps
# @btime from measuring global-variable access.
let p = ahmadfard_cases(4)
    Q = Float64.(p.Q); xy = p.xy; V = p.V / p.nb
    rb, D, ks, Cs, s = p.r.b, p.D, p.k.s, p.C.s, p.s
    ro, ri, kg, kp = p.r.o, p.r.i, p.k.g, p.k.p
    kf, cf, ρf, μf = p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ
    T0, Tlim = p.T.g, [p.T.L, p.T.H]
    println("\nSizing time per level (case 4):")
    for lvl in levels
        print("  ", lvl, ": ")
        @btime alternative_sizing($Q, $xy, $rb, $D, $ks, $Cs, $s, $ro, $ri, $kg, $kp, $kf,
            $cf, $ρf, $μf, $V, $T0, $Tlim; level = $lvl)
    end
end

# Computed vs reference lengths, grouped by case.
fig = Figure(size = (760, 460))
ax = Axis(fig[1, 1], xlabel = "Test case", ylabel = "Borehole length H [m]",
    title = "Alternative ASHRAE sizing — computed vs Dion & Pasquier (2025) reference",
    xticks = (collect(cases), string.(cases)))
colors = Makie.wong_colors()
for (li, lvl) in enumerate(levels)
    x = collect(cases) .+ (li - 2) * 0.22
    barplot!(ax, x, H[:, li]; width = 0.22, color = colors[li], label = string(lvl))
    scatter!(ax, x, H_ref[:, li]; color = :black, marker = :hline, markersize = 18)
end
axislegend(ax, "Level (bars); ─ = reference"; position = :lt)
display(fig)
