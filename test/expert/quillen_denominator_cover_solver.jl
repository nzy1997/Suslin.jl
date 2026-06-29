using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_PATCH_CATALOG_PATH)
end

function _solver_tamper_multiplier(result)
    multipliers = copy(result.coverage_multipliers)
    multipliers[1] += one(result.ring)
    return Suslin.QuillenDenominatorCoverSolverResult(
        result.source_candidate,
        result.ring,
        result.raw_denominators,
        result.exponent,
        result.powered_denominators,
        multipliers,
        result.coverage_terms,
        result.coverage_sum,
        result.cover_certificate,
        result.verification,
    )
end

function _solver_tamper_exponent(result)
    return Suslin.QuillenDenominatorCoverSolverResult(
        result.source_candidate,
        result.ring,
        result.raw_denominators,
        result.exponent + 1,
        result.powered_denominators,
        result.coverage_multipliers,
        result.coverage_terms,
        result.coverage_sum,
        result.cover_certificate,
        result.verification,
    )
end

function _coverage_error_message(err)
    return err isa ArgumentError && occursin("coverage not proven", sprint(showerror, err))
end

@testset "Quillen denominator cover solver" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()

    two_open = entries["quillen-two-open-cover-qq"]
    R = two_open.ring.object
    raw = [data.denominator for data in two_open.denominator_data]
    result = Suslin.solve_quillen_denominator_cover(R, raw; max_exponent = 2)

    @test result isa Suslin.QuillenDenominatorCoverSolverResult
    @test Suslin.verify_quillen_denominator_cover_solver_result(result)
    @test Suslin.verify_quillen_denominator_cover(result.cover_certificate)
    @test result.raw_denominators == raw
    @test result.exponent == 1
    @test result.powered_denominators == raw
    @test result.coverage_multipliers == [one(R), one(R)]
    @test result.coverage_terms == [result.coverage_multipliers[i] * result.powered_denominators[i] for i in eachindex(raw)]
    @test result.coverage_sum == one(R)
    @test result.cover_certificate.denominators == result.powered_denominators
    replay = Suslin.replay_quillen_denominator_cover_solver_result(result)
    @test replay.powered_denominators == result.powered_denominators
    @test replay.coverage_sum == result.coverage_sum
    @test replay.coverage_ok == result.verification.coverage_ok

    squared = Suslin.solve_quillen_denominator_cover(
        R,
        raw;
        max_exponent = 2,
        exponent = 2,
    )
    @test squared.exponent == 2
    @test squared.powered_denominators == [denominator^2 for denominator in raw]
    @test squared.raw_denominators == raw
    @test squared.cover_certificate.denominators == squared.powered_denominators
    @test Suslin.verify_quillen_denominator_cover_solver_result(squared)

    nontrivial = entries["quillen-nontrivial-multipliers-qq"]
    R_nontrivial = nontrivial.ring.object
    raw_nontrivial = [data.denominator for data in nontrivial.denominator_data]
    supplied = [data.coverage_multiplier for data in nontrivial.denominator_data]
    supplied_result = Suslin.solve_quillen_denominator_cover(
        R_nontrivial,
        raw_nontrivial;
        max_exponent = 1,
        coverage_multipliers = supplied,
    )
    @test supplied_result.coverage_multipliers == supplied
    @test any(multiplier -> multiplier != one(R_nontrivial), supplied_result.coverage_multipliers)
    @test supplied_result.coverage_sum == one(R_nontrivial)
    @test Suslin.verify_quillen_denominator_cover_solver_result(supplied_result)

    @test !Suslin.verify_quillen_denominator_cover_solver_result(_solver_tamper_multiplier(result))
    @test !Suslin.verify_quillen_denominator_cover_solver_result(_solver_tamper_exponent(result))

    R_bad, (X_bad, r_bad, s_bad) = Oscar.polynomial_ring(QQ, ["X", "r", "s"])
    try
        Suslin.solve_quillen_denominator_cover(R_bad, [r_bad, s_bad]; max_exponent = 2)
        @test false
    catch err
        @test _coverage_error_message(err)
    end

    try
        Suslin.solve_quillen_denominator_cover(
            R_nontrivial,
            raw_nontrivial;
            max_exponent = 1,
            coverage_multipliers = [supplied[1] + one(R_nontrivial), supplied[2]],
        )
        @test false
    catch err
        @test _coverage_error_message(err)
    end
end
