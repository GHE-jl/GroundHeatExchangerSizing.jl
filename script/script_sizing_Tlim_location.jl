# Impact of where the temperature limits Tlim are applied in the alternative ASHRAE sizing. The
# locations (keyword Tlim_at) are seen from the ground, not the heat pump:
#   :inlet (default) — fluid entering the ground -> heat-pump leaving fluid temperature (LFT)
#   :mean — mean fluid temperature in the borehole
#   :outlet — fluid leaving the ground -> heat-pump entering fluid temperature (EFT)
# The four Ahmadfard & Bernier (2019) cases are sized at the three levels (L2, L3, L4) for each
# location, and the lengths are printed in a table.

using Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate()

using GroundHeatExchangerSizing

# Reference cases (hourly ground loads + all parameters) from Ahmadfard & Bernier (2019).
include("Ahmadfard_cases.jl")

levels = (:L2, :L3, :L4)
locations = (:inlet, :mean, :outlet)
cases = 1:4

H = zeros(length(cases), length(levels), length(locations))
for (ci, c) in enumerate(cases)
    p = ahmadfard_cases(c)
    Q = Float64.(p.Q)
    V = p.V / p.nb                     # per-borehole loop flow [m³/s] (case stores total flow)
    for (li, lvl) in enumerate(levels), (ai, at) in enumerate(locations)
        H[ci, li, ai] = alternative_sizing(Q, p.xy, p.r.b, p.D, p.k.s, p.C.s, p.s, p.r.o, p.r.i,
            p.k.g, p.k.p, p.k.f, p.C.f / p.ρ.f, p.ρ.f, p.μ, V, p.T.g, [p.T.L, p.T.H];
            level = lvl, Tlim_at = at)
    end
end

# Borehole lengths [m] and the difference of :inlet and :outlet relative to :mean [%]
println("Case  Level   inlet (m)   mean (m)   outlet (m)   inlet vs mean (%)   outlet vs mean (%)")
for (ci, c) in enumerate(cases), (li, lvl) in enumerate(levels)
    Hi, Hm, Ho = H[ci, li, :]
    println(lpad(c, 4), lpad(string(lvl), 7), lpad(round(Hi; digits = 1), 12),
        lpad(round(Hm; digits = 1), 11), lpad(round(Ho; digits = 1), 13),
        lpad(round(100 * (Hi - Hm) / Hm; digits = 1), 20),
        lpad(round(100 * (Ho - Hm) / Hm; digits = 1), 21))
end
