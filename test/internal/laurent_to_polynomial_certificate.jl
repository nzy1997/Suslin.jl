using Test
using Suslin
using Oscar

function _ltp_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _ltp_apply_factors(factors, column, R)
    return _ltp_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _ltp_matrix_column_to_tuple(column_matrix)
    @test ncols(column_matrix) == 1
    return tuple((column_matrix[row, 1] for row in 1:nrows(column_matrix))...)
end

function _ltp_tamper(certificate; kwargs...)
    fields = propertynames(certificate)
    values = map(field -> get(kwargs, field, getproperty(certificate, field)), fields)
    return Suslin.LaurentToPolynomialColumnCertificate(values...)
end

function _ltp_tamper_noether(certificate; kwargs...)
    fields = propertynames(certificate)
    values = map(field -> get(kwargs, field, getproperty(certificate, field)), fields)
    return Suslin.LaurentNoetherCertificate(values...)
end

function _ltp_captured_error(callback)
    try
        callback()
        return nothing
    catch err
        return err
    end
end

function _ltp_modify_first_factor_coefficient(certificate)
    R = certificate.ring
    n = length(certificate.original_column)
    modified = collect(certificate.conversion_factors)
    idx = findfirst(factor -> Suslin._canonical_elementary_factor_record(factor).kind == :elementary, modified)
    idx === nothing && error("expected an elementary factor to modify")
    record = Suslin._canonical_elementary_factor_record(modified[idx])
    modified[idx] = elementary_matrix(n, record.row, record.col, record.coefficient + one(R), R)
    return modified
end

function _ltp_assert_common_certificate(certificate)
    @test certificate.validation_status == :ok
    @test Suslin._validate_laurent_to_polynomial_certificate(certificate) == :ok
    @test certificate.replay.noether_status == :ok
    @test certificate.replay.elementary_factors_are_elementary
    @test certificate.replay.elementary_replay_ok
    @test certificate.replay.polynomial_entries_in_ring
    @test certificate.replay.inverse_lift_ok
    @test certificate.replay.polynomial_unimodular
    @test certificate.replay.overall_ok

    R = certificate.ring
    P = certificate.polynomial_ring
    n = length(certificate.original_column)
    @test length(certificate.original_column) >= 3
    @test certificate.conversion_product == _ltp_factor_product(certificate.conversion_factors, R, n)
    replayed = _ltp_matrix_column_to_tuple(_ltp_apply_factors(
        certificate.conversion_factors,
        certificate.elementary_source_column,
        R,
    ))
    @test replayed == certificate.intermediate_laurent_column
    @test all(factor -> Suslin._is_elementary_matrix_factor(factor, R, n), certificate.conversion_factors)
    @test all(entry -> parent(entry) === P, certificate.polynomial_column)
    @test Suslin.is_unimodular_column(collect(certificate.polynomial_column), P)
    @test tuple((
        Suslin._laurent_to_polynomial_lift_entry(entry, certificate.inverse_lift)
        for entry in certificate.polynomial_column
    )...) == certificate.intermediate_laurent_column

    ordinary_factor = elementary_matrix(n, 1, 2, gens(P)[1] + gens(P)[2], P)
    lifted_factor = Suslin._laurent_to_polynomial_lift_factor(
        ordinary_factor,
        certificate.factor_lift_metadata,
    )
    @test base_ring(lifted_factor) === R
    @test lifted_factor[1, 2] == gens(R)[1] + gens(R)[2]
    @test all(row == col || lifted_factor[row, col] == zero(R) || (row == 1 && col == 2)
        for row in 1:n, col in 1:n)
end

@testset "Laurent-to-polynomial column conversion certificate" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    polynomial_column = [x, y, one(R) + x]
    polynomial_noether = Suslin._laurent_noether_certificate(polynomial_column, 1, x)
    polynomial_certificate = Suslin._laurent_to_polynomial_certificate(
        polynomial_column,
        polynomial_noether,
        1,
        x,
    )
    _ltp_assert_common_certificate(polynomial_certificate)
    @test isempty(polynomial_certificate.conversion_factors)
    @test polynomial_certificate.elementary_source_column == polynomial_noether.transformed_column
    @test polynomial_certificate.intermediate_laurent_column == polynomial_noether.transformed_column
    @test tuple(string.(gens(polynomial_certificate.polynomial_ring))...) == ("x", "y")

    U, (u, v) = suslin_laurent_polynomial_ring(QQ, ["u", "v"])
    unit_column = [u^-1, v^-1 + u, one(U) + v]
    unit_noether = Suslin._laurent_noether_certificate(unit_column, 1, u)
    unit_certificate = Suslin._laurent_to_polynomial_certificate(unit_column, unit_noether, 1, u)
    _ltp_assert_common_certificate(unit_certificate)
    @test length(unit_certificate.conversion_factors) >= 2
    @test unit_certificate.intermediate_laurent_column == (zero(U), zero(U), one(U))
    @test tuple(string.(unit_certificate.polynomial_column)...) == ("0", "0", "1")
    @test any(entry -> any(exp -> any(<(0), exp), collect(exponents(entry))), unit_certificate.elementary_source_column)

    reordered_factors = reverse(unit_certificate.conversion_factors)
    @test Suslin._validate_laurent_to_polynomial_certificate(_ltp_tamper(
        unit_certificate;
        conversion_factors = reordered_factors,
        conversion_product = _ltp_factor_product(reordered_factors, U, length(unit_column)),
    )) != :ok

    modified_factors = _ltp_modify_first_factor_coefficient(unit_certificate)
    @test Suslin._validate_laurent_to_polynomial_certificate(_ltp_tamper(
        unit_certificate;
        conversion_factors = modified_factors,
        conversion_product = _ltp_factor_product(modified_factors, U, length(unit_column)),
    )) != :ok

    @test Suslin._validate_laurent_to_polynomial_certificate(_ltp_tamper(
        unit_certificate;
        intermediate_laurent_column = (one(U), zero(U), zero(U)),
    )) != :ok

    bad_forward = (
        merge(unit_certificate.forward_polynomialization[1], (; polynomial_generator = unit_certificate.polynomial_generators[2])),
        unit_certificate.forward_polynomialization[2],
    )
    @test Suslin._validate_laurent_to_polynomial_certificate(_ltp_tamper(
        unit_certificate;
        forward_polynomialization = bad_forward,
    )) != :ok

    bad_inverse = (
        merge(unit_certificate.inverse_lift[1], (; laurent_value = gens(U)[2])),
        unit_certificate.inverse_lift[2],
    )
    @test Suslin._validate_laurent_to_polynomial_certificate(_ltp_tamper(
        unit_certificate;
        inverse_lift = bad_inverse,
    )) != :ok

    P = unit_certificate.polynomial_ring
    @test Suslin._validate_laurent_to_polynomial_certificate(_ltp_tamper(
        unit_certificate;
        polynomial_column = (one(P), zero(P), zero(P)),
    )) != :ok

    stale_noether = _ltp_tamper_noether(
        unit_noether;
        transformed_column = tuple(unit_column...),
    )
    @test Suslin._validate_laurent_to_polynomial_certificate(_ltp_tamper(
        unit_certificate;
        noether_certificate = stale_noether,
    )) != :ok
    @test_throws ArgumentError Suslin._laurent_to_polynomial_certificate(
        unit_column,
        stale_noether,
        1,
        u,
    )

    nonunimodular_column = [u + one(U), (u + one(U)) * v, zero(U)]
    nonunimodular_noether = Suslin._laurent_noether_certificate(nonunimodular_column, 1, u)
    @test !Suslin.is_unimodular_column(nonunimodular_column, U)
    @test_throws ArgumentError Suslin._laurent_to_polynomial_certificate(
        nonunimodular_column,
        nonunimodular_noether,
        1,
        u,
    )

    unsupported_column = [u^-1 + one(U), v + one(U), u + v]
    unsupported_noether = Suslin._laurent_noether_certificate(unsupported_column, 1, u)
    @test Suslin.is_unimodular_column(unsupported_column, U)
    @test Suslin._validate_laurent_noether_certificate(unsupported_noether) == :ok
    @test !any(is_unit, unsupported_noether.transformed_column)
    @test any(entry -> any(exp -> any(<(0), Int.(collect(exp))), collect(exponents(entry))),
        unsupported_noether.transformed_column)
    unsupported_error = _ltp_captured_error(() -> Suslin._laurent_to_polynomial_certificate(
        unsupported_column,
        unsupported_noether,
        1,
        u,
    ))
    @test unsupported_error isa ArgumentError
    @test occursin("must already be polynomial or contain a Laurent unit entry", sprint(showerror, unsupported_error))
end
