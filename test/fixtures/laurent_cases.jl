module LaurentFixtureCatalog

using Oscar
using Suslin

include("toricbuilder_factor_toric_block_3.jl")

function _column_matrix(R, values)
    M = zero_matrix(R, length(values), 1)
    for (i, value) in enumerate(values)
        M[i, 1] = value
    end
    return M
end

function _ring_metadata(R, x, y)
    return (;
        description = "GF(2)[x^+/-1, y^+/-1]",
        object = R,
        generators = (x, y),
    )
end

function _ring_constructor_metadata()
    return (;
        function_name = :suslin_laurent_polynomial_ring,
        coefficient = "GF(2)",
        variables = ("x", "y"),
    )
end

function _synthetic_provenance(description)
    return (;
        source = :synthetic,
        issue = "#6",
        description,
    )
end

function _pinv_toricbuilder_case(ring, ring_constructor)
    toricbuilder_fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    pinv = only(filter(entry -> entry.toricbuilder_role == "Pinv", toricbuilder_fixture.cases))

    return (;
        id = "toricbuilder-factor-toric-block-3-pinv",
        kind = :toricbuilder_relation,
        ring_constructor,
        ring,
        dimensions = (;
            matrix = pinv.size,
            source_matrix = pinv.size,
        ),
        inputs = (;
            matrix = pinv.matrix,
            source_matrix = pinv.source_matrix,
        ),
        expected_relation = (;
            kind = :toricbuilder_left_inverse,
            description = pinv.relation_description,
        ),
        provenance = (;
            source = :toricbuilder_contract_fixture,
            issue = "#19",
            fixture_id = pinv.name,
            toricbuilder_role = pinv.toricbuilder_role,
            toricbuilder_commit = pinv.provenance.toricbuilder_commit,
            generation_command = pinv.provenance.generation_command,
        ),
        determinant_profile = (;
            relevant = true,
            expected_class = pinv.determinant_classification,
        ),
        consumer_test_ids = ("issue-6-laurent-fixtures", "issue-19-toricbuilder-contract"),
    )
end

function catalog()
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    ring = _ring_metadata(R, x, y)
    ring_constructor = _ring_constructor_metadata()

    solvable_matrix = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])
    solvable_solution = _column_matrix(R, [one(R), y^-1])
    solvable_rhs = solvable_matrix * solvable_solution

    unsolvable_matrix = zero_matrix(R, 1, 1)
    unsolvable_rhs = _column_matrix(R, [one(R)])

    negative_vector = _column_matrix(R, [x^-1 * y, x^-2])
    normalization_unit = x^2
    normalized_vector = _column_matrix(R, [x * y, one(R)])

    return (;
        ring,
        cases = [
            (;
                id = "laurent-linear-system-solvable",
                kind = :solvable_linear_system,
                ring_constructor,
                ring,
                dimensions = (;
                    matrix = (2, 2),
                    rhs = (2, 1),
                    solution = (2, 1),
                ),
                inputs = (;
                    matrix = solvable_matrix,
                    rhs = solvable_rhs,
                    expected_solution = solvable_solution,
                ),
                expected_relation = (;
                    kind = :linear_system_solution,
                    description = "matrix * expected_solution == rhs",
                ),
                provenance = _synthetic_provenance("triangular 2x2 Laurent linear system with exact solution"),
                determinant_profile = (;
                    relevant = true,
                    expected_class = "one",
                ),
                consumer_test_ids = ("issue-6-laurent-fixtures", "issue-8-linear-laurent-tests"),
            ),
            (;
                id = "laurent-linear-system-unsolvable",
                kind = :unsolvable_linear_system,
                ring_constructor,
                ring,
                dimensions = (;
                    matrix = (1, 1),
                    rhs = (1, 1),
                ),
                inputs = (;
                    matrix = unsolvable_matrix,
                    rhs = unsolvable_rhs,
                    unsolvability_certificate = (;
                        kind = :zero_matrix_nonzero_rhs,
                        rhs_index = (1, 1),
                    ),
                ),
                expected_relation = (;
                    kind = :linear_system_no_solution,
                    description = "zero matrix cannot produce a nonzero right-hand side",
                ),
                provenance = _synthetic_provenance("zero 1x1 Laurent linear system with nonzero rhs"),
                determinant_profile = (;
                    relevant = true,
                    expected_class = "non-unit",
                ),
                consumer_test_ids = ("issue-6-laurent-fixtures", "issue-9-linear-laurent-rejections"),
            ),
            (;
                id = "laurent-negative-exponent-normalization",
                kind = :negative_exponent_normalization,
                ring_constructor,
                ring,
                dimensions = (;
                    vector = (2, 1),
                    normalized_vector = (2, 1),
                ),
                inputs = (;
                    vector = negative_vector,
                    normalization_unit,
                    normalized_vector,
                ),
                expected_relation = (;
                    kind = :negative_exponent_normalization,
                    description = "x^2 * [x^-1*y, x^-2] == [x*y, 1]",
                ),
                provenance = _synthetic_provenance("minimal vector normalization with negative x exponents"),
                determinant_profile = (;
                    relevant = false,
                    expected_class = nothing,
                ),
                consumer_test_ids = ("issue-6-laurent-fixtures", "issue-12-negative-exponent-normalization"),
            ),
            _pinv_toricbuilder_case(ring, ring_constructor),
        ],
    )
end

end
