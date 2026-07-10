using Test
using Suslin
using Oscar

const TORICBUILDER_CACHE_CASE010_REPORT_SCRIPT =
    joinpath(@__DIR__, "..", "..", "scripts", "report_toricbuilder_cache_q_blocks.jl")

include(TORICBUILDER_CACHE_CASE010_REPORT_SCRIPT)

function _case010_entry()
    return only(filter(
        entry -> entry.id == "case_010",
        ToricBuilderCacheQBlockStatusReport.ToricBuilderCacheQBlocks.catalog().cases,
    ))
end

function _case010_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case010_certificate_with_original(certificate, original_matrix)
    return Suslin.LaurentGLFactorizationCertificate(
        original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

function _case010_corrupt_matrix(A)
    corrupted = copy(A)
    corrupted[1, 1] += one(base_ring(A))
    return corrupted
end

function _case010_corrupt_entry(entry)
    corrupted_entries = collect(entry.sparse_entries)
    row, col, _ = corrupted_entries[1]
    corrupted_entries[1] = (row, col, "0")
    return merge(entry, (; sparse_entries = corrupted_entries, sparse_entry_count = length(corrupted_entries)))
end

function _case010_last_column_at_dimension(A, target_dimension::Int)
    current = normalize_laurent_gl_matrix(A).normalized_matrix
    while nrows(current) > target_dimension
        current = Suslin._laurent_column_peel_step(current).next_block
    end
    nrows(current) == target_dimension || error("case_010 dimension $(target_dimension) was not reached")
    return [current[row, target_dimension] for row in 1:target_dimension]
end

function _case010_verify_or_staged_false(certificate)
    try
        return verify_laurent_gl_factorization_certificate(certificate)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        message = sprint(showerror, err)
        occursin("staged", message) || occursin("unsupported", message) || rethrow()
        return false
    end
end

@testset "ToricBuilder case_010 Laurent GL certificate" begin
    entry = _case010_entry()
    A = ToricBuilderCacheQBlockStatusReport.ToricBuilderCacheQBlocks.materialize_matrix(entry)
    R = base_ring(A)

    profile = classify_laurent_determinant(A)
    @test profile.classification == :laurent_monomial_unit
    laurent_gens = collect(gens(R))
    @test profile.determinant == laurent_gens[1] * laurent_gens[2]

    certificate = laurent_gl_factorization_certificate(A)
    @test certificate.original_matrix == A
    @test certificate.normalized_core == certificate.normalization.normalized_matrix
    @test det(certificate.normalized_core) == one(R)
    @test length(certificate.core_factors) > 0
    @test length(certificate.core_factorization.peel_steps) == nrows(A) - 2
    @test any(step -> step.dimension == 4, certificate.core_factorization.peel_steps)
    @test _case010_product(certificate.core_factors, R, nrows(A)) == certificate.normalized_core
    @test certificate.reconstructed_product == A
    @test verify_laurent_gl_factorization_certificate(certificate)

    dimension4_column = _case010_last_column_at_dimension(A, 4)
    dimension4_candidate = Suslin._laurent_unit_creation_candidate(dimension4_column, R)
    @test dimension4_candidate.pivot_index == 3
    @test dimension4_candidate.source_index == 1

    public_error = try
        elementary_factorization(A)
        nothing
    catch err
        err
    end
    @test public_error isa ArgumentError
    public_message = sprint(showerror, public_error)
    @test occursin("elementary_factorization(A) is an elementary-only SL_n API", public_message)
    @test occursin("laurent_gl_factorization_certificate(A)", public_message)

    corrupted_matrix = _case010_corrupt_matrix(A)
    corrupted_certificate = _case010_certificate_with_original(certificate, corrupted_matrix)
    @test !_case010_verify_or_staged_false(corrupted_certificate)

    corrupted_entry = _case010_corrupt_entry(entry)
    corrupted_row = ToricBuilderCacheQBlockStatusReport._exercised_row(corrupted_entry)
    @test !(corrupted_row.route_status == :gl_certificate_pass && corrupted_row.verified)
end
