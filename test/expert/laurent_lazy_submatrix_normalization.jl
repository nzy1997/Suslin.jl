using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue157_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue157_deferred_probe(certificate, records)
    return function (candidate)
        push!(records, (; size = (nrows(candidate), ncols(candidate)), candidate))
        @test candidate == certificate.deferred_submatrix
        @test candidate != certificate.original_matrix
        return Suslin.classify_laurent_determinant(candidate)
    end
end

function _issue157_embedded_deferred(deferred)
    R = base_ring(deferred)
    return block_embedding(deferred, nrows(deferred) + 1, collect(1:nrows(deferred)))
end

@testset "deferred Laurent submatrix determinant-one normalization metadata" begin
    entry = _issue157_fixture("determinant-one-triangular")
    certificate = Suslin._laurent_determinant_deferred_peel_certificate(entry.inputs.matrix)
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(certificate)
    R = base_ring(certificate.deferred_submatrix)

    @test metadata.peel_certificate == certificate
    @test metadata.determinant_source == :deferred_submatrix
    @test metadata.overall_determinant == one(R)
    @test metadata.determinant_classification == :one
    @test metadata.supported
    @test metadata.deferred_correction.kind == :identity
    @test metadata.deferred_diagonal_correction === nothing
    @test metadata.normalized_deferred_core == certificate.deferred_submatrix
    @test det(metadata.normalized_deferred_core) == one(R)
    @test metadata.staged_boundary === nothing
    @test metadata.verification.overall_ok
    @test Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)
end

@testset "deferred Laurent submatrix monomial-unit normalization metadata" begin
    R, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    deferred = matrix(R, [
        u * v one(R);
        zero(R) one(R)
    ])
    A = _issue157_embedded_deferred(deferred)
    certificate = Suslin._laurent_determinant_deferred_peel_certificate(A)
    probe_records = Any[]
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(
        certificate;
        determinant_probe = _issue157_deferred_probe(certificate, probe_records),
    )

    @test length(probe_records) == 1
    @test only(probe_records).size == (2, 2)
    @test metadata.overall_determinant == u * v
    @test metadata.determinant_classification == :laurent_monomial_unit
    @test metadata.supported
    @test metadata.deferred_correction.kind == :left_diagonal_determinant_correction
    @test metadata.deferred_diagonal_correction == metadata.deferred_correction
    @test metadata.deferred_correction.factor * metadata.normalized_deferred_core ==
        certificate.deferred_submatrix
    @test det(metadata.normalized_deferred_core) == one(R)
    @test metadata.staged_boundary === nothing
    @test metadata.verification.overall_ok
    @test metadata.verification.normalized_core_det_ok
    @test Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)
end

@testset "deferred Laurent submatrix non-unit staged boundary" begin
    entry = _issue157_fixture("non-unit-determinant-negative")
    A = _issue157_embedded_deferred(entry.inputs.matrix)
    certificate = Suslin._laurent_determinant_deferred_peel_certificate(A)
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(certificate)

    @test certificate.deferred_submatrix == entry.inputs.matrix
    @test metadata.overall_determinant == entry.determinant_profile.expected_determinant
    @test metadata.determinant_classification == :non_unit
    @test !metadata.supported
    @test metadata.deferred_correction === nothing
    @test metadata.deferred_diagonal_correction === nothing
    @test metadata.normalized_deferred_core === nothing
    @test metadata.staged_boundary !== nothing
    @test metadata.staged_boundary.kind == :unsupported_deferred_laurent_determinant
    @test metadata.staged_boundary.reason == :non_unit_deferred_determinant
    @test metadata.staged_boundary.overall_determinant == metadata.overall_determinant
    @test metadata.verification.overall_ok
    @test metadata.verification.boundary_ok
    @test Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)
end

@testset "deferred Laurent submatrix unsupported unit staged boundary" begin
    Z4, _ = residue_ring(ZZ, 4)
    S, (u,) = laurent_polynomial_ring(Z4, ["u"])
    deferred = matrix(S, [
        one(S) + S(Z4(2)) * u zero(S);
        zero(S) one(S)
    ])
    @test classify_laurent_determinant(deferred).classification == :other_unit

    A = _issue157_embedded_deferred(deferred)
    certificate = Suslin._laurent_determinant_deferred_peel_certificate(A)
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(certificate)

    @test metadata.determinant_classification == :other_unit
    @test !metadata.supported
    @test metadata.deferred_correction === nothing
    @test metadata.deferred_diagonal_correction === nothing
    @test metadata.normalized_deferred_core === nothing
    @test metadata.staged_boundary.kind == :unsupported_deferred_laurent_determinant
    @test metadata.staged_boundary.reason == :unsupported_deferred_unit_class
    @test metadata.verification.boundary_ok
    @test Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)
end

@testset "lazy GL peel returns enriched metadata for non-one deferred determinant" begin
    entry = _issue157_fixture("monomial-unit-row-column-cores")
    metadata = Suslin._factor_laurent_gl_lazy_determinant_peel(entry.inputs.matrix)

    @test metadata.determinant_source == :deferred_submatrix
    @test metadata.overall_determinant == entry.determinant_profile.expected_determinant
    @test metadata.determinant_classification == :laurent_monomial_unit
    @test metadata.supported
    @test metadata.normalized_deferred_core !== nothing
    @test det(metadata.normalized_deferred_core) == one(base_ring(metadata.normalized_deferred_core))
    @test metadata.staged_boundary === nothing
end
