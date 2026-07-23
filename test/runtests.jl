using GroundHeatExchangerSizing
using Test

# Synthetic, perfectly balanced hourly ground load (single borehole)
# A daily+seasonal sine with zero annual mean, so both operating limits are exercised.
const NH = 8760
const _t_hours = collect(1:NH)
const Q_synthetic = 4000.0 .* sinpi.(2 .* _t_hours ./ NH) .+ 800.0 .* sinpi.(2 .* _t_hours ./ 24)

# Single-borehole geometry and properties (loosely Case 1 of Ahmadfard & Bernier 2019).
const XY1 = reshape([0.0, 0.0], 1, 2)
# Cg and Cp set to the DeepANN-valid values it would otherwise silently clamp them to (Cg to its
# [1.4e6, 3.0e6] upper bound, Cp to its fixed 1.9e6) — see `outlet_sizing`'s `DeepANN` default.
const CASE = (rb = 0.075, D = 4.0, ks = 1.8, Cs = 2.0736e6, s = 0.075, ro = 0.0167, ri = 0.0137,
    kg = 1.4, Cg = 3.0e6, kp = 0.43, Cp = 1.9e6, kf = 0.48, cf = 3795.0, ρf = 1052.0,
    μf = 5.2e-3, V = 4.0e-4, T0 = 17.5, Tlim = [0.0, 35.0])

# Helper: assert an outlet sizing result is well-formed and inside the search bounds.
function check_outlet(H)
    @test H isa Real
    @test isfinite(H)
    @test GroundHeatExchangerSizing._H_LB - 1e-6 <= H <= GroundHeatExchangerSizing._H_UB + 1e-6
end

# Helper: assert an alternative sizing result is well-formed.
function check_alternative(H)
    @test H isa Real
    @test isfinite(H)
end

@testset "GroundHeatExchangerSizing.jl" begin

    @testset "Thermal load analysis" begin
        Qm = Q_hourly_to_monthly(Q_synthetic)
        @test size(Qm) == (12, 3)
        @test all(Qm[:, 2] .<= 0.0)        # monthly cooling peaks are non-positive
        @test all(Qm[:, 3] .>= 0.0)        # monthly heating peaks are non-negative

        Q3 = Q_hourly_to_three_pulses(Q_synthetic)
        @test size(Q3) == (3, 2)
        @test Q3[1, 1] == Q3[1, 2]         # yearly average shared by both columns
        @test Q3[3, 1] ≈ minimum(Q_synthetic)
        @test Q3[3, 2] ≈ maximum(Q_synthetic)

        Qhr = Q_monthly_to_hourly(Qm, 6 * 3600.0)
        @test size(Qhr) == (NH, 2)

        Q3m = Q_monthly_to_three_pulses(Qm)
        @test size(Q3m) == (3, 2)

        analysis = Q_analysis(Q_synthetic)
        @test analysis isa QLoads
        @test analysis.h ≈ maximum(abs.(Q_synthetic))

        # Q_cutoff scales peaks; use a copy since it mutates its argument.
        Qc = Q_cutoff(copy(Q_synthetic), 0.5, 0.5)
        @test maximum(Qc) ≈ 0.5 * maximum(Q_synthetic)
    end

    @testset "Alternative ASHRAE sizing" begin
        # Reference (Ahmadfard & Bernier 2019, Appendix B; identical inputs): validation deferred.
        for level in (:L2, :L3, :L4)
            res = alternative_sizing(Q_synthetic, XY1, CASE.rb, CASE.D, CASE.ks, CASE.Cs, CASE.s,
                CASE.ro, CASE.ri, CASE.kg, CASE.kp, CASE.kf, CASE.cf, CASE.ρf, CASE.μf, CASE.V,
                CASE.T0, CASE.Tlim; level = level, ny = 2.0)
            check_alternative(res)
        end
        # Explicit level entry points agree with the dispatcher.
        r_disp = alternative_sizing(Q_synthetic, XY1, CASE.rb, CASE.D, CASE.ks, CASE.Cs, CASE.s,
            CASE.ro, CASE.ri, CASE.kg, CASE.kp, CASE.kf, CASE.cf, CASE.ρf, CASE.μf, CASE.V,
            CASE.T0, CASE.Tlim; level = :L2, ny = 2.0)
        r_L2 = alternative_sizing_L2(Q_hourly_to_three_pulses(Q_synthetic), XY1, CASE.rb, CASE.D,
            CASE.ks, CASE.Cs, CASE.s, CASE.ro, CASE.ri, CASE.kg, CASE.kp, CASE.kf, CASE.cf,
            CASE.ρf, CASE.μf, CASE.V, CASE.T0, CASE.Tlim; ny = 2.0)
        @test r_disp ≈ r_L2

        @test_throws ArgumentError alternative_sizing(Q_synthetic, XY1, CASE.rb, CASE.D, CASE.ks,
            CASE.Cs, CASE.s, CASE.ro, CASE.ri, CASE.kg, CASE.kp, CASE.kf, CASE.cf, CASE.ρf,
            CASE.μf, CASE.V, CASE.T0, CASE.Tlim; level = :L5)
    end

    @testset "Borehole-outlet transfer-function sizing" begin
        # Quantitative comparison against Dion & Pasquier 2025, Table 2 is done on the paper's own
        # four cases in script/script_outlet_sizing.jl; this synthetic load only checks the result
        # is well-formed.
        for level in (:L2, :L3, :L4)
            res = outlet_sizing(Q_synthetic, XY1, CASE.rb, CASE.D, CASE.ks, CASE.Cs, CASE.s,
                CASE.ro, CASE.ri, CASE.kg, CASE.Cg, CASE.kp, CASE.Cp, CASE.kf, CASE.cf, CASE.ρf,
                CASE.μf, CASE.V, CASE.T0, CASE.Tlim; level = level, ny = 2.0)
            check_outlet(res)
        end

        @test_throws ArgumentError outlet_sizing(Q_synthetic, XY1, CASE.rb, CASE.D, CASE.ks,
            CASE.Cs, CASE.s, CASE.ro, CASE.ri, CASE.kg, CASE.Cg, CASE.kp, CASE.Cp, CASE.kf,
            CASE.cf, CASE.ρf, CASE.μf, CASE.V, CASE.T0, CASE.Tlim; level = :L5)
    end
end
