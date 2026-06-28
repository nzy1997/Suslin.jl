using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue159_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue159_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue159_rebuild(
    certificate;
    determinant_source = certificate.determinant_source,
    correction_side = certificate.correction_side,
    reconstruction_relation = certificate.reconstruction_relation,
    correction = certificate.correction,
    inverse_correction = certificate.inverse_correction,
    normalized_deferred_core = certificate.normalized_deferred_core,
    normalized_deferred_factorization = certificate.normalized_deferred_factorization,
    normalized_deferred_factors = certificate.normalized_deferred_factors,
    rewritten_left_factors = certificate.rewritten_left_factors,
    rewritten_right_factors = certificate.rewritten_right_factors,
    elementary_factors = certificate.elementary_factors,
    elementary_product = certificate.elementary_product,
    reconstructed_product = certificate.reconstructed_product,
    verification = certificate.verification,
)
    return Suslin.LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
        certificate.overall_determinant,
        determinant_source,
        correction_side,
        reconstruction_relation,
        correction,
        inverse_correction,
        normalized_deferred_core,
        normalized_deferred_factorization,
        normalized_deferred_factors,
        rewritten_left_factors,
        rewritten_right_factors,
        elementary_factors,
        elementary_product,
        reconstructed_product,
        verification,
    )
end

function _issue159_caught_error(f)
    try
        f()
        return nothing
    catch caught
        return caught
    end
end

@testset "lazy Laurent determinant correction side choices" begin
    entry = _issue159_fixture("issue-38-q-block-lazy-determinant")
    A = entry.inputs.matrix
    R = base_ring(A)
    n = nrows(A)

    row_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :row,
    )
    column_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :column,
    )

    @test row_certificate.correction_side == :row
    @test column_certificate.correction_side == :column
    @test row_certificate.correction.side == :left
    @test column_certificate.correction.side == :right
    @test row_certificate.reconstruction_relation == :left_correction_times_elementary_product
    @test column_certificate.reconstruction_relation == :elementary_product_times_right_correction
    @test row_certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test column_certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test row_certificate.overall_determinant == column_certificate.overall_determinant

    @test row_certificate.elementary_product ==
        _issue159_product(row_certificate.elementary_factors, R, n)
    @test column_certificate.elementary_product ==
        _issue159_product(column_certificate.elementary_factors, R, n)
    @test row_certificate.correction.factor * row_certificate.elementary_product == A
    @test column_certificate.elementary_product * column_certificate.correction.factor == A
    @test row_certificate.reconstructed_product == A
    @test column_certificate.reconstructed_product == A
    @test row_certificate.verification.overall_ok
    @test column_certificate.verification.overall_ok
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(row_certificate)
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(column_certificate)

    @test det(row_certificate.correction.factor) == row_certificate.overall_determinant
    @test det(column_certificate.correction.factor) == column_certificate.overall_determinant
    @test row_certificate.correction.factor * row_certificate.inverse_correction ==
        identity_matrix(R, n)
    @test column_certificate.correction.factor * column_certificate.inverse_correction ==
        identity_matrix(R, n)

    wrong_metadata = _issue159_rebuild(
        column_certificate;
        correction_side = :row,
        reconstruction_relation = :left_correction_times_elementary_product,
    )
    wrong_verification = Suslin._laurent_gl_lazy_deferred_correction_certificate_verification(
        wrong_metadata,
    )
    @test !wrong_verification.overall_ok
    @test !wrong_verification.correction_side_ok
    @test !wrong_verification.correction_ok
    @test !wrong_verification.rewritten_left_factors_ok
    @test !Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(wrong_metadata)

    err = try
        Suslin._laurent_gl_lazy_deferred_correction_certificate(
            A;
            correction_side = :diagonal,
        )
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin(":row", sprint(showerror, err))
    @test occursin(":column", sprint(showerror, err))
end

@testset "lazy Laurent determinant-one correction side choices" begin
    entry = _issue159_fixture("determinant-one-triangular")
    A = entry.inputs.matrix
    R = base_ring(A)
    n = nrows(A)

    row_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :row,
    )
    column_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :column,
    )

    @test row_certificate.overall_determinant == one(R)
    @test column_certificate.overall_determinant == one(R)
    @test row_certificate.correction.kind == :identity
    @test column_certificate.correction.kind == :identity
    @test row_certificate.correction.side == :left
    @test column_certificate.correction.side == :right
    @test row_certificate.normalized_deferred_core ==
        row_certificate.deferred_metadata.deferred_submatrix
    @test column_certificate.normalized_deferred_core ==
        column_certificate.deferred_metadata.deferred_submatrix
    @test row_certificate.correction.factor * row_certificate.elementary_product == A
    @test column_certificate.elementary_product * column_certificate.correction.factor == A
    @test row_certificate.correction.factor == identity_matrix(R, n)
    @test column_certificate.correction.factor == identity_matrix(R, n)
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(row_certificate)
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(column_certificate)
end

@testset "lazy Laurent correction side guard coverage" begin
    entry = _issue159_fixture("monomial-unit-row-column-cores")
    metadata = Suslin._factor_laurent_gl_lazy_determinant_peel(entry.inputs.matrix)
    certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        metadata;
        correction_side = :row,
    )

    external_err = _issue159_caught_error(() ->
        Suslin._laurent_gl_external_correction_side(:diagonal))
    @test external_err isa ArgumentError
    @test occursin("unsupported internal Laurent correction side", sprint(showerror, external_err))

    relation_err = _issue159_caught_error(() ->
        Suslin._laurent_gl_reconstruction_relation(:diagonal))
    @test relation_err isa ArgumentError
    @test occursin("unsupported internal Laurent correction side", sprint(showerror, relation_err))

    unsupported_metadata = merge(metadata, (; determinant_classification = :other_unit))
    unsupported_err = _issue159_caught_error(() ->
        Suslin._laurent_lazy_deferred_correction(unsupported_metadata, :left))
    @test unsupported_err isa ArgumentError
    @test occursin("unsupported deferred Laurent determinant", sprint(showerror, unsupported_err))

    bad_side_correction = merge(metadata.deferred_correction, (; side = :middle))
    normalized_err = _issue159_caught_error(() ->
        Suslin._laurent_lazy_normalized_deferred_core(metadata, bad_side_correction))
    @test normalized_err isa ArgumentError
    @test occursin("unsupported internal Laurent correction side", sprint(showerror, normalized_err))

    bad_embedded_correction = merge(certificate.correction, (; side = :middle))
    assembly_err = _issue159_caught_error(() ->
        Suslin._laurent_lazy_hoist_elementary_factors(
            metadata,
            bad_embedded_correction,
            certificate.normalized_deferred_factors,
        ))
    @test assembly_err isa ArgumentError
    @test occursin("unsupported internal Laurent correction side", sprint(showerror, assembly_err))
end
