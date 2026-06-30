using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_PATCH_CATALOG_PATH)
end

function _solver_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function _solver_local_correction(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function _solver_sequence_certificate(entry, index::Int)
    local_factor = entry.local_factors[index]
    realization = Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = _solver_local_certificate(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = _solver_local_correction(local_factor),
        factors = [local_factor.factor],
        patched_substitution_witness = entry.patched_substitution_witness,
        ring = entry.ring.object,
        size = entry.size,
    )
    return Suslin.quillen_local_factor_sequence_certificate(
        realization;
        factor_provenance = (;
            factor_index = 1,
            sequence_index = 1,
            local_index = 1,
            fixture_id = entry.id,
        ),
        metadata = (;
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )
end

function _solver_denominator_cover_candidate(entry)
    return Suslin.extract_quillen_denominator_cover_candidate([
        _solver_sequence_certificate(entry, index)
        for index in eachindex(entry.local_factors)
    ])
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

    candidate = _solver_denominator_cover_candidate(two_open)
    candidate_result = Suslin.solve_quillen_denominator_cover(candidate; max_exponent = 2)
    @test candidate_result.source_candidate === candidate
    @test candidate_result.raw_denominators == candidate.raw_denominators
    @test candidate_result.powered_denominators == candidate.raw_denominators
    @test candidate_result.verification.source_candidate_ok
    @test Suslin.verify_quillen_denominator_cover_solver_result(candidate_result)

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
    malformed_solver_result = Suslin.QuillenDenominatorCoverSolverResult(
        result.source_candidate,
        QQ,
        result.raw_denominators,
        result.exponent,
        result.powered_denominators,
        result.coverage_multipliers,
        result.coverage_terms,
        result.coverage_sum,
        result.cover_certificate,
        result.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover_solver_result(malformed_solver_result)

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
