module ToricBuilderFactorToricBlock3Fixture

using Oscar

function _sparse_laurent_matrix(R, rows::Int, cols::Int, entries)
    M = zero_matrix(R, rows, cols)
    for (i, j, value) in entries
        M[i, j] = value
    end
    return M
end

function fixture()
    R, (x, y) = Oscar.laurent_polynomial_ring(Oscar.GF(2), ["x", "y"])

    qinv = _sparse_laurent_matrix(R, 16, 16, [
        (1, 1, one(R)), (1, 3, x), (1, 5, one(R)), (1, 6, y),
        (2, 2, one(R)), (2, 4, x), (2, 5, one(R)), (2, 6, one(R)),
        (3, 3, x + 1), (3, 5, one(R)), (3, 6, y), (3, 7, one(R)), (3, 8, y),
        (4, 4, one(R)), (5, 6, one(R)), (6, 4, one(R)), (6, 7, one(R)),
        (7, 3, one(R)), (7, 4, one(R)), (8, 6, one(R)), (8, 8, one(R)),
        (9, 9, one(R)), (10, 10, one(R)), (11, 9, one(R)), (11, 10, one(R)), (11, 13, one(R)),
        (12, 11, one(R)), (12, 12, one(R)), (12, 13, x^-1), (12, 15, one(R)),
        (13, 9, y^-1), (13, 10, one(R)), (13, 14, one(R)), (13, 16, one(R)),
        (14, 9, one(R)), (14, 10, one(R)), (14, 13, one(R)), (14, 15, one(R)),
        (15, 9, one(R)), (15, 10, 1 + x^-1), (15, 11, one(R)), (15, 13, 1 + x^-1),
        (16, 9, y^-1), (16, 10, y^-1), (16, 13, y^-1), (16, 16, one(R)),
    ])

    column_transformation = _sparse_laurent_matrix(R, 16, 16, [
        (1, 1, one(R)), (1, 3, one(R)), (1, 5, y), (1, 6, one(R)), (1, 7, one(R)), (1, 8, y),
        (2, 2, one(R)), (2, 3, one(R)), (2, 5, one(R)), (2, 6, one(R)), (2, 7, x + 1), (2, 8, y),
        (3, 4, one(R)), (3, 7, one(R)), (4, 4, one(R)),
        (5, 3, one(R)), (5, 4, x), (5, 6, one(R)), (5, 7, x + 1), (5, 8, y),
        (6, 5, one(R)), (7, 4, one(R)), (7, 6, one(R)), (8, 5, one(R)), (8, 8, one(R)),
        (9, 9, one(R)), (10, 10, one(R)), (11, 9, x^-1), (11, 11, 1 + x^-1), (11, 15, one(R)),
        (12, 10, x^-1), (12, 12, one(R)), (12, 14, one(R)), (12, 15, one(R)),
        (13, 9, one(R)), (13, 10, one(R)), (13, 11, one(R)),
        (14, 9, y^-1), (14, 10, one(R)), (14, 11, y^-1), (14, 13, one(R)), (14, 16, one(R)),
        (15, 11, one(R)), (15, 14, one(R)), (16, 11, y^-1), (16, 16, one(R)),
    ])

    pinv = _sparse_laurent_matrix(R, 8, 8, [
        (1, 1, one(R)), (2, 2, one(R)), (3, 1, one(R)), (3, 3, one(R)),
        (4, 2, one(R)), (4, 3, one(R)), (4, 4, one(R)),
        (5, 7, one(R)), (6, 6, one(R)), (7, 5, one(R)),
        (8, 5, one(R)), (8, 6, x^-1), (8, 7, one(R)), (8, 8, one(R)),
    ])

    row_transformation = _sparse_laurent_matrix(R, 8, 8, [
        (1, 1, one(R)), (2, 2, one(R)), (3, 1, one(R)), (3, 3, one(R)),
        (4, 1, one(R)), (4, 2, one(R)), (4, 3, one(R)), (4, 4, one(R)),
        (5, 7, one(R)), (6, 6, one(R)), (7, 5, one(R)),
        (8, 5, one(R)), (8, 6, x^-1), (8, 7, one(R)), (8, 8, one(R)),
    ])

    ring = (;
        description = "GF(2)[x^+/-1, y^+/-1]",
        base_field = "GF(2)",
        variables = ("x", "y"),
        kind = "multivariate Laurent polynomial ring",
    )
    provenance = (;
        toricbuilder_path = "/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder",
        toricbuilder_commit = "fa7f82252d42fdc0b2726bc48af24ac4c70a8d73",
        source_function = "src/toric_form/toric_factorization.jl:factor_toric_block",
        generation_command = "factor_toric_block(3, x, y, R) with R, (x, y) = laurent_polynomial_ring(GF(2), [\"x\", \"y\"])",
    )

    return (;
        ring,
        block_size = (2, 2),
        provenance,
        cases = [
            (;
                name = "factor_toric_block_3_qinv",
                toricbuilder_role = "Qinv",
                matrix = qinv,
                source_matrix = column_transformation,
                relation_description = "column_transformation * Qinv == I_16",
                ring = ring.description,
                size = (16, 16),
                determinant_classification = "one",
                expected_suslin_status = :supported_column_peel,
                expected_suslin_path = :laurent_column_peel,
                expected_output = :verified_transformation_certificate,
                provenance,
            ),
            (;
                name = "factor_toric_block_3_pinv",
                toricbuilder_role = "Pinv",
                matrix = pinv,
                source_matrix = row_transformation,
                relation_description = "row_transformation * Pinv == I_8",
                ring = ring.description,
                size = (8, 8),
                determinant_classification = "one",
                expected_suslin_status = :supported_column_peel,
                expected_suslin_path = :laurent_column_peel,
                expected_output = :verified_transformation_certificate,
                provenance,
            ),
        ],
    )
end

end
