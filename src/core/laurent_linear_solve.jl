function _validate_laurent_linear_inputs(A, B)
    R = _require_laurent_polynomial_ring(base_ring(A); label="A base ring")

    base_ring(B) === R || throw(ArgumentError("Matrices A and B must be over the same Laurent polynomial ring"))

    rows_A = nrows(A)
    rows_B = nrows(B)
    rows_B == rows_A || throw(DimensionMismatch("Number of rows in A ($(rows_A)) must match number of rows in B ($(rows_B))"))

    return R
end

function _solve_laurent_linear_native(A, B)
    solvable, solution = can_solve_with_solution(A, B; side=:right)
    solvable && return solution
    throw(ErrorException("No exact solution exists for A * U = B"))
end

function _is_native_laurent_solver_capability_error(err, backtrace)
    frames = stacktrace(backtrace)
    touches_native_solver = any(frame -> startswith(String(frame.func), "can_solve_with_solution"), frames)
    touches_native_solver || return false

    if err isa ErrorException
        return sprint(showerror, err) == "divrem not implemented for LaurentMPoly"
    end

    err isa MethodError || return false

    missing_name = string(err.f)
    (missing_name == "gcdxx" || missing_name == "annihilator") || return false

    return true
end

function _laurent_linear_module_data(A, quotient_map)
    A_quo = map_entries(quotient_map, A)
    quotient_ring = base_ring(A_quo)
    polynomial_ring = base_ring(quotient_ring)
    quotient_ideal = modulus(quotient_ring)
    A_poly = map_entries(x -> polynomial_ring(Oscar.lift(x)), A_quo)

    free = free_module(polynomial_ring, nrows(A_poly))
    coefficient_count = ncols(A_poly)
    generators = [free(collect(A_poly[:, j])) for j in 1:coefficient_count]

    row_count = nrows(A_poly)
    for generator in gens(quotient_ideal)
        for row in 1:row_count
            relation = [zero(polynomial_ring) for _ in 1:row_count]
            relation[row] = generator
            push!(generators, free(relation))
        end
    end

    submodule, _ = sub(free, generators)
    return (;
        free,
        submodule,
        coefficient_count,
        polynomial_ring,
        quotient_ring,
    )
end

function _solve_laurent_linear_column(b_col, module_data, quotient_map, R)
    (; free, submodule, coefficient_count, polynomial_ring, quotient_ring) = module_data

    b_quo = map_entries(quotient_map, b_col)
    b_poly = map_entries(x -> polynomial_ring(Oscar.lift(x)), b_quo)
    b_vec = free(collect(b_poly[:, 1]))
    coords = Oscar.coordinates(b_vec, submodule)

    entries = [
        preimage(quotient_map, quotient_ring(coords[i]))
        for i in 1:coefficient_count
    ]
    return matrix(R, length(entries), 1, entries)
end

function _solve_laurent_linear_fallback(A, B)
    R = _validate_laurent_linear_inputs(A, B)
    quotient_map = Oscar._polyringquo(R)
    rhs_cols = ncols(B)
    module_data = _laurent_linear_module_data(A, quotient_map)

    try
        rhs_cols == 1 && return _solve_laurent_linear_column(B, module_data, quotient_map, R)

        solutions = [
            _solve_laurent_linear_column(B[:, j:j], module_data, quotient_map, R)
            for j in 1:rhs_cols
        ]
        return hcat(solutions...)
    catch err
        if err isa ErrorException && occursin("not liftable to the given generating system", sprint(showerror, err))
            throw(ErrorException("No exact solution exists for A * U = B"))
        end
        rethrow()
    end
end

"""
    solve_laurent_linear(A, B)

Solve the Laurent-polynomial linear system `A * U = B` exactly.

The implementation prefers Oscar's native right-side linear solver. If the
installed Oscar/AbstractAlgebra stack does not support Laurent-polynomial
matrices for that path, it falls back to a module-based exact solver while
preserving the same public API and no-solution error contract.

# Arguments
- `A`: `m x n` matrix over a Laurent polynomial ring
- `B`: `m x k` matrix over the same Laurent polynomial ring

# Returns
- `U`: `n x k` solution matrix such that `A * U == B`

# Throws
- `ErrorException` if no exact solution exists
- `ArgumentError` if the inputs are not over the same Laurent polynomial ring
- `DimensionMismatch` if `nrows(B) != nrows(A)`
"""
function solve_laurent_linear(A, B)
    _validate_laurent_linear_inputs(A, B)

    try
        return _solve_laurent_linear_native(A, B)
    catch err
        if _is_native_laurent_solver_capability_error(err, catch_backtrace())
            return _solve_laurent_linear_fallback(A, B)
        end
        rethrow()
    end
end
