using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")

function _issue160_caught_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "public Laurent GL certificate options" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    Q = entry.inputs.matrix
    expected_determinant = det(Q)

    eager_certificate = laurent_gl_factorization_certificate(Q)
    @test eager_certificate isa LaurentGLFactorizationCertificate
    @test eager_certificate.reconstructed_product == Q
    @test verify_laurent_gl_factorization_certificate(eager_certificate)

    explicit_eager = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :eager,
    )
    @test explicit_eager isa LaurentGLFactorizationCertificate
    @test explicit_eager.reconstructed_product == Q
    @test verify_laurent_gl_factorization_certificate(explicit_eager)

    row_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    @test row_certificate isa LaurentLazyGLHoistCertificate
    @test row_certificate.overall_determinant == expected_determinant
    @test row_certificate.determinant_source == :deferred_submatrix
    @test row_certificate.correction_side == :row
    @test row_certificate.reconstructed_product == Q
    @test row_certificate.verification.overall_ok
    @test verify_laurent_gl_factorization_certificate(row_certificate)

    default_row_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
    )
    @test default_row_certificate.correction_side == :row
    @test default_row_certificate.determinant_source == :deferred_submatrix
    @test verify_laurent_gl_factorization_certificate(default_row_certificate)

    column_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :column,
    )
    @test column_certificate isa LaurentLazyGLHoistCertificate
    @test column_certificate.overall_determinant == expected_determinant
    @test column_certificate.overall_determinant == row_certificate.overall_determinant
    @test column_certificate.determinant_source == :deferred_submatrix
    @test column_certificate.correction_side == :column
    @test column_certificate.reconstructed_product == Q
    @test column_certificate.verification.overall_ok
    @test verify_laurent_gl_factorization_certificate(column_certificate)

    strategy_err = _issue160_caught_error(() ->
        laurent_gl_factorization_certificate(Q; determinant_strategy = :deferred))
    @test strategy_err isa ArgumentError
    @test occursin(":eager", sprint(showerror, strategy_err))
    @test occursin(":lazy", sprint(showerror, strategy_err))

    misplaced_side_err = _issue160_caught_error(() ->
        laurent_gl_factorization_certificate(Q; correction_side = :column))
    @test misplaced_side_err isa ArgumentError
    @test occursin("determinant_strategy = :lazy", sprint(showerror, misplaced_side_err))

    invalid_side_err = _issue160_caught_error(() ->
        laurent_gl_factorization_certificate(
            Q;
            determinant_strategy = :lazy,
            correction_side = :diagonal,
        ))
    @test invalid_side_err isa ArgumentError
    @test occursin(":row", sprint(showerror, invalid_side_err))
    @test occursin(":column", sprint(showerror, invalid_side_err))

    original_err = _issue160_caught_error(() -> elementary_factorization(Q))
    @test original_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, original_err))
end
