using Test
using Suslin
using Oscar

function exact_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function exact_apply_factors(factors, v, R)
    n = length(v)
    column = matrix(R, n, 1, collect(v))
    return exact_reduction_product(factors, R, n) * column
end

function exact_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function assert_reduces_to_last_unit(v, R)
    factors = Suslin.reduce_unimodular_column(v, R)
    @test exact_apply_factors(factors, v, R) == exact_target_column(R, length(v))
    return factors
end

function assert_factors_reduce_to_last_unit(factors, v, R)
    @test factors !== nothing
    @test exact_apply_factors(factors, v, R) == exact_target_column(R, length(v))
    return factors
end

function captured_reduction_error(v, R)
    try
        Suslin.reduce_unimodular_column(v, R)
        return nothing
    catch err
        return err
    end
end

@testset "exact unimodular reduction supports longer ordinary columns" begin
    F2 = GF(2)
    R, (x, y) = Oscar.polynomial_ring(F2, ["x", "y"])

    hard_slice = [
        x + y^2,
        x * y + x + one(R),
        x^2 + x * y + y + one(R),
    ]

    length6 = vcat(hard_slice, [x^2, x * y, y^2 + x])
    @test Suslin.is_unimodular_column(length6, R)
    assert_reduces_to_last_unit(length6, R)

    length8 = [
        x^2 + y,
        hard_slice[1],
        x * y,
        y^2 + one(R),
        hard_slice[2],
        x^3 + y,
        x * y + y,
        hard_slice[3],
    ]
    @test Suslin.is_unimodular_column(length8, R)
    assert_reduces_to_last_unit(length8, R)
    assert_factors_reduce_to_last_unit(Suslin._reduce_via_supported_three_block(length8, R), length8, R)

    witness_slice = [x, y, x + one(R)]
    length12 = [
        x^2,
        y^2,
        x * y,
        witness_slice[1],
        x^2 + y,
        y^3 + x,
        witness_slice[2],
        x * y + y,
        x^3 + y^2,
        y^2 + y,
        x^2 * y + x,
        witness_slice[3],
    ]
    @test Suslin.is_unimodular_column(length12, R)
    assert_reduces_to_last_unit(length12, R)
end

@testset "exact unimodular reduction supports Laurent-normalized columns" begin
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])

    v = [
        x^-1,
        x^-2 * y,
        x^-1 + x^-2,
        x^-1 * y,
        x^-2 * y^2 + x^-1,
        x^-2,
    ]

    @test Suslin.is_unimodular_column(v, R)
    assert_reduces_to_last_unit(v, R)

    normalized = [
        x^-1 + x^-2 * y^2,
        x^-1 * y + x^-1 + x^-2,
        one(R) + x^-1 * y + x^-2 * y + x^-2,
        x^-1 + x^-2 * y,
        x^-1 * y + x^-2 * y^2,
        x^-2 * y + x^-1 * y^2,
    ]
    @test !any(is_unit, normalized)
    @test Suslin.is_unimodular_column(normalized, R)
    assert_reduces_to_last_unit(normalized, R)

    non_unimodular_normalized = [x + y, x * y + y, x^2 + x * y]
    @test !Suslin.is_unimodular_column(non_unimodular_normalized, R)
    @test Suslin._reduce_laurent_unimodular_column(non_unimodular_normalized, R) === nothing

    unsupported_normalized = [x * y + x, x^2 + x + one(R), x * y + y^2 + one(R)]
    @test Suslin.is_unimodular_column(unsupported_normalized, R)
    @test Suslin._reduce_laurent_unimodular_column(unsupported_normalized, R) === nothing
end

@testset "exact unimodular reduction preserves old small cases" begin
    F2 = GF(2)
    R2, (x, y) = Oscar.polynomial_ring(F2, ["x", "y"])

    assert_reduces_to_last_unit([x, y, one(R2)], R2)
    assert_reduces_to_last_unit([x, one(R2), y], R2)
    assert_reduces_to_last_unit([x, y, x + one(R2)], R2)
    assert_reduces_to_last_unit([
        x + y^2,
        x * y + x + one(R2),
        x^2 + x * y + y + one(R2),
    ], R2)
    monicity_fallback = [
        x^2 + x * y^2 + x + y + one(R2),
        x^2 + x * y^2 + x * y + x + y^2 + one(R2),
        x^2 * y + x * y + one(R2),
    ]
    @test Suslin._reduce_after_monicity_normalization_certificate(monicity_fallback, R2) !== nothing
    @test_throws ArgumentError Suslin._reduce_via_general_ecp_pipeline_certificate(monicity_fallback, R2)
    fallback_cert = Suslin.ecp_column_reduction_certificate(monicity_fallback, R2)
    @test Suslin.verify_ecp_column_reduction(fallback_cert)
    @test any(stage -> stage.kind == :monicity_normalization, fallback_cert.stages)
    fallback_diagnostic = Suslin.diagnose_unimodular_column_reduction(monicity_fallback, R2)
    @test fallback_diagnostic.status == :supported
    @test :general_ecp_pipeline in fallback_diagnostic.attempted_stages
    @test :monicity_normalization in fallback_diagnostic.attempted_stages

    F3 = GF(3)
    R3, (t,) = Oscar.polynomial_ring(F3, ["t"])
    assert_reduces_to_last_unit([t + one(R3), t, R3(2)], R3)
    assert_reduces_to_last_unit([t + one(R3), R3(2), t], R3)
end

@testset "exact unimodular reduction staged failures" begin
    F2 = GF(2)
    R, (x, y) = Oscar.polynomial_ring(F2, ["x", "y"])

    non_unimodular = [x, y, x * y]
    non_unimodular_err = captured_reduction_error(non_unimodular, R)
    @test non_unimodular_err isa ArgumentError
    @test occursin("v must be a unimodular column", sprint(showerror, non_unimodular_err))

    unsupported = [zero(R), x^2, x * y + one(R)]
    @test Suslin.is_unimodular_column(unsupported, R)
    unsupported_err = captured_reduction_error(unsupported, R)
    @test unsupported_err isa ArgumentError
    @test occursin("unsupported exact unimodular column reduction", sprint(showerror, unsupported_err))
    @test occursin("general ECP pipeline", sprint(showerror, unsupported_err))
    @test !occursin("not unimodular", sprint(showerror, unsupported_err))

    short_err = captured_reduction_error([one(R), zero(R)], R)
    @test short_err isa ArgumentError
    @test occursin("length at least 3", sprint(showerror, short_err))

    @test !Suslin._has_at_least_two_generators(:not_a_ring)
    @test Suslin._reduce_via_supported_three_block([x, y, x + one(R)], R) === nothing
    @test Suslin._reduce_via_supported_three_block([x, x * y, y, x^2 * y], R) === nothing
end
