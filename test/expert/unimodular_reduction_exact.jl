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
    @test !occursin("not unimodular", sprint(showerror, unsupported_err))
end
