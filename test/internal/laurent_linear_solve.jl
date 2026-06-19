using Test
using Suslin
using Oscar

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
