using Test
using Suslin
using Oscar

const ISSUE9_STUB_RING = first(suslin_laurent_polynomial_ring(GF(2), ["issue9_stub"]))

function _issue9_column_matrix(R, values)
    M = zero_matrix(R, length(values), 1)
    for (i, value) in enumerate(values)
        M[i, 1] = value
    end
    return M
end

function _issue9_error_message(err)
    return sprint(showerror, err)
end

struct _Issue9NativeRethrowMatrix
    rows::Int
    cols::Int
end

struct _Issue9FallbackRethrowMatrix
    rows::Int
    cols::Int
    ring
end

AbstractAlgebra.base_ring(::_Issue9NativeRethrowMatrix) = ISSUE9_STUB_RING
AbstractAlgebra.nrows(A::_Issue9NativeRethrowMatrix) = A.rows
AbstractAlgebra.ncols(A::_Issue9NativeRethrowMatrix) = A.cols
AbstractAlgebra.can_solve_with_solution(::_Issue9NativeRethrowMatrix, ::_Issue9NativeRethrowMatrix; side=:right) =
    throw(ErrorException("unexpected native solver failure"))

AbstractAlgebra.base_ring(A::_Issue9FallbackRethrowMatrix) = A.ring
AbstractAlgebra.nrows(A::_Issue9FallbackRethrowMatrix) = A.rows
AbstractAlgebra.ncols(A::_Issue9FallbackRethrowMatrix) = A.cols
Suslin._solve_laurent_linear_column(::_Issue9FallbackRethrowMatrix, module_data, quotient_map, R) =
    throw(ErrorException("unexpected fallback failure"))

@testset "Laurent linear solve" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])

    A = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])

    single_solution = _issue9_column_matrix(R, [one(R), y^-1])
    single_rhs = A * single_solution
    computed_single = solve_laurent_linear(A, single_rhs)
    @test computed_single == single_solution
    @test A * computed_single == single_rhs

    multi_solution = matrix(R, [
        one(R)      y;
        x^-1        one(R) + y^-1
    ])
    multi_rhs = A * multi_solution
    computed_multi = solve_laurent_linear(A, multi_rhs)
    @test computed_multi == multi_solution
    @test A * computed_multi == multi_rhs

    unsolvable_A = zero_matrix(R, 1, 1)
    unsolvable_B = _issue9_column_matrix(R, [one(R)])
    try
        solve_laurent_linear(unsolvable_A, unsolvable_B)
        error("expected no-solution failure")
    catch err
        @test err isa ErrorException
        @test _issue9_error_message(err) == "No exact solution exists for A * U = B"
    end

    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    wrong_parent_rhs = _issue9_column_matrix(S, [one(S), v])
    @test_throws ArgumentError solve_laurent_linear(A, wrong_parent_rhs)

    wrong_rows_rhs = _issue9_column_matrix(R, [one(R)])
    @test_throws DimensionMismatch solve_laurent_linear(A, wrong_rows_rhs)
end

@testset "Laurent linear solve internal error boundaries" begin
    P, _ = suslin_polynomial_ring(QQ, ["t"])
    native_A = matrix(P, 1, 1, [one(P)])
    native_B = matrix(P, 1, 1, [one(P)])
    @test Suslin._solve_laurent_linear_native(native_A, native_B) == native_B

    try
        Suslin._solve_laurent_linear_native(zero_matrix(P, 1, 1), native_B)
        error("expected native no-solution failure")
    catch err
        @test err isa ErrorException
        @test _issue9_error_message(err) == "No exact solution exists for A * U = B"
    end

    native_failure = _Issue9NativeRethrowMatrix(1, 1)
    try
        solve_laurent_linear(native_failure, native_failure)
        error("expected native rethrow")
    catch err
        @test err isa ErrorException
        @test _issue9_error_message(err) == "unexpected native solver failure"
    end

    R, (x, _) = suslin_laurent_polynomial_ring(GF(2), ["fallback_x", "fallback_y"])
    fallback_A = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])
    fallback_rhs = _Issue9FallbackRethrowMatrix(2, 1, R)
    try
        Suslin._solve_laurent_linear_fallback(fallback_A, fallback_rhs)
        error("expected fallback rethrow")
    catch err
        @test err isa ErrorException
        @test _issue9_error_message(err) == "unexpected fallback failure"
    end
end
