using Oscar
using Suslin

const ISSUE38_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "..", "test", "fixtures", "toricbuilder_issue38_cases.jl")

if !isdefined(Main, :ToricBuilderIssue38Cases)
    include(ISSUE38_FIXTURE_PATH)
end

function _factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _report_original_boundary(Q)
    err = try
        factors = elementary_factorization(Q)
        error("original Issue 38 Q unexpectedly factored with $(length(factors)) factors")
    catch caught
        caught
    end

    err isa ArgumentError || throw(err)
    message = sprint(showerror, err)
    occursin("Laurent GL_n normalization boundary", message) ||
        error("unexpected original Q factorization error: $(message)")

    println("original_q status=STAGED_BOUNDARY det=$(det(Q))")
    println("original_q error=\"$(message)\"")
    return err
end

function _run_row_unit_correction(Q)
    normalization = normalize_laurent_gl_matrix(Q)
    core = normalization.normalized_matrix
    R = base_ring(core)
    n = nrows(core)
    factors = elementary_factorization(core)
    product = _factor_product(factors, R, n)

    product == core || error("row correction factors do not multiply to the normalized core")
    verify_factorization(core, factors) || error("row correction factorization did not verify")
    normalization.correction.factor * product == Q ||
        error("row correction factor and core product do not reconstruct Q")

    println(
        "row_unit_correction status=PASS det_core=$(det(core)) factors=$(length(factors)) verified=true",
    )
    return (; core, factors, correction = normalization.correction.factor)
end

function _run_column_unit_correction(Q)
    R = base_ring(Q)
    n = nrows(Q)
    determinant = det(Q)
    correction = diagonal_matrix(R, [inv(determinant), [one(R) for _ in 2:n]...])
    inverse_correction = diagonal_matrix(R, [determinant, [one(R) for _ in 2:n]...])
    core = Q * correction
    factors = elementary_factorization(core)
    product = _factor_product(factors, R, n)

    det(core) == one(R) || error("column correction core is not determinant one")
    product == core || error("column correction factors do not multiply to the normalized core")
    verify_factorization(core, factors) || error("column correction factorization did not verify")
    product * inverse_correction == Q ||
        error("column correction core product and inverse correction do not reconstruct Q")

    println(
        "column_unit_correction status=PASS det_core=$(det(core)) factors=$(length(factors)) verified=true",
    )
    return (; core, factors, correction, inverse_correction)
end

function _run_lazy_gl_certificate(Q)
    certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    verified = verify_laurent_gl_factorization_certificate(certificate)

    certificate.reconstructed_product == Q ||
        error("lazy GL certificate does not reconstruct original Q")
    verified || error("lazy GL certificate did not verify")

    println(
        "lazy_gl_certificate status=PASS det=$(certificate.overall_determinant) correction_side=$(certificate.correction_side) determinant_source=$(certificate.determinant_source) factors=$(length(certificate.elementary_factors)) verified=$(verified)",
    )
    return certificate
end

function main()
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    Q = entry.inputs.matrix
    profile = classify_laurent_determinant(Q)

    println("issue38_q size=$(size(Q)) determinant=$(profile.determinant) class=$(profile.classification)")
    _report_original_boundary(Q)
    row = _run_row_unit_correction(Q)
    column = _run_column_unit_correction(Q)
    lazy = _run_lazy_gl_certificate(Q)
    println(
        "issue38_unit_correction_example status=PASS row_factors=$(length(row.factors)) column_factors=$(length(column.factors)) lazy_factors=$(length(lazy.elementary_factors))",
    )
    return (; Q, row, column, lazy)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
