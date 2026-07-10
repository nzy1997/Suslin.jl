module LaurentFixtureCatalog

using Oscar
using Suslin

include("toricbuilder_factor_toric_block_3.jl")
include("toricbuilder_case008_d14_column_boundary.jl")

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

function _laurent_to_poly_source_refs()
    return (;
        laurent_to_poly = (;
            author = "Park",
            algorithm = :algorithm_6_1,
            name = "LaurentToPoly",
            role = :laurent_to_polynomial_conversion,
        ),
        laurent_noether = (;
            author = "Park",
            algorithm = :algorithm_6_3,
            name = "LaurentNoether",
            role = :laurent_variable_change_normalization,
        ),
    )
end

function _laurent_to_poly_route_ring_metadata(R, variables)
    return (;
        description = "GF(2)[$(variables[1])^+/-1, $(variables[2])^+/-1]",
        object = R,
        generators = variables,
    )
end

function _laurent_to_poly_route_constructor(variables)
    return (;
        function_name = :suslin_laurent_polynomial_ring,
        coefficient = "GF(2)",
        variables = Tuple(string.(variables)),
    )
end

function _laurent_to_poly_entry_term_count(entry)
    iszero(entry) && return 0
    return length(collect(coefficients(entry)))
end

function _laurent_to_poly_selected_entry_contract(
    source_column,
    selected_entry_index,
    selected_entry_role,
)
    selected_entry = source_column[selected_entry_index]
    return (;
        preserves_unimodularity = true,
        polynomial_target = true,
        selected_entry_must_be_polynomial_unit = true,
        selected_entry_index,
        selected_entry_role,
        selected_source_fingerprint =
            Suslin._laurent_descent_column_support_fingerprint([selected_entry]),
        selected_source_term_count = _laurent_to_poly_entry_term_count(selected_entry),
        selected_source_is_unit = is_unit(selected_entry),
    )
end

function _laurent_to_poly_route_entry(
    id,
    route,
    R,
    variables,
    source_column,
    selected_entry_index,
    expected_reducer,
    provenance_extra,
    verifier_id,
    selected_entry_role,
)
    return (;
        id,
        kind = :laurent_to_polynomial_route,
        route,
        ring_constructor = _laurent_to_poly_route_constructor(variables),
        ring = _laurent_to_poly_route_ring_metadata(R, variables),
        source_column,
        selected_entry_index,
        source_fingerprint = Suslin._laurent_descent_column_support_fingerprint(source_column),
        expected_reducer,
        provenance = (;
            route,
            source_refs = _laurent_to_poly_source_refs(),
            provenance_extra...,
        ),
        post_conversion_contract = _laurent_to_poly_selected_entry_contract(
            source_column,
            selected_entry_index,
            selected_entry_role,
        ),
        verifier_id,
        consumer_test_ids = (
            "issue-351-laurent-to-poly-route-fixtures",
            "issue-351-laurent-to-poly-route-consumer",
        ),
    )
end

function laurent_to_poly_route_catalog()
    normalization_ring, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    normalization_column = [
        x^-1 * y^-1 * (x + y^2),
        x^-1 * y^-1 * (x * y + x + one(normalization_ring)),
        x^-1 * y^-1 * (x^2 + x * y + y + one(normalization_ring)),
        zero(normalization_ring),
        zero(normalization_ring),
        zero(normalization_ring),
    ]

    general_ring, (general_x, general_y) =
        Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    general_column = [
        general_x * general_y + general_x,
        general_x^2 + general_x + one(general_ring),
        general_x * general_y + general_y^2 + one(general_ring),
    ]

    d14 = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(d14) == :ok ||
        throw(ArgumentError("case008 d14 boundary fixture failed validation"))
    return (;
        cases = [
            _laurent_to_poly_route_entry(
                "laurent-to-poly-existing-normalization",
                :existing_normalization,
                normalization_ring,
                (x, y),
                normalization_column,
                1,
                (;
                    status = :supported,
                    failure_code = nothing,
                    stage = :laurent_normalization,
                    stage_outcome = :delegated_to_polynomial,
                    normalized_status = :supported,
                    normalized_failure_code = nothing,
                ),
                (; source_case = :laurent_column_reduction_diagnostics),
                :laurent_to_poly_existing_normalization,
                :existing_normalization_anchor,
            ),
            _laurent_to_poly_route_entry(
                "laurent-to-poly-general-ecp",
                :general_ecp,
                general_ring,
                (general_x, general_y),
                general_column,
                2,
                (;
                    status = :unsupported,
                    failure_code = :unsupported_laurent_column_family,
                    stage = :laurent_native_ecp_boundary,
                    stage_outcome = :staged_boundary,
                    boundary = :laurent_native_ecp,
                    requires_descent_measure = true,
                    certified_descent_scope = nothing,
                    next_boundary = nothing,
                    requires_link_witness = true,
                    requires_endpoint_reduction = true,
                    requires_laurent_normality_replay = true,
                    requires_recursive_peel_integration = true,
                    fallback_policy = :diagnostic_only,
                ),
                (; source_case = :laurent_column_reduction_diagnostics),
                :laurent_to_poly_general_ecp,
                :general_ecp_anchor,
            ),
            _laurent_to_poly_route_entry(
                "laurent-to-poly-case008-d14",
                :case008_d14,
                d14.ring,
                Tuple(gens(d14.ring)),
                d14.failing_column,
                1,
                (;
                    status = :unsupported,
                    failure_code = :unsupported_laurent_column_family,
                    stage = :laurent_native_ecp_boundary,
                    stage_outcome = :staged_boundary,
                    boundary = :laurent_native_ecp,
                    requires_descent_measure = false,
                    certified_descent_scope = :single_certified_step,
                    next_boundary = :laurent_endpoint_reduction,
                    requires_link_witness = false,
                    requires_endpoint_reduction = true,
                    requires_laurent_normality_replay = true,
                    requires_recursive_peel_integration = true,
                    fallback_policy = :diagnostic_only,
                ),
                (;
                    source_case = d14.case_id,
                    source_fixture = :ToricBuilderCase008D14ColumnBoundary,
                    case_id = d14.case_id,
                    dimension = d14.first_failing_peel_dimension,
                    source_cache_file = d14.source_cache_file,
                    source_block = d14.source_block,
                    source_matrix_dimensions = d14.source_matrix_dimensions,
                    source_column_transformation_dimensions =
                        d14.source_column_transformation_dimensions,
                    passed_peel_dimensions = d14.passed_peel_dimensions,
                    first_failing_peel_dimension = d14.first_failing_peel_dimension,
                    boundary_status = :current_staged_d14_boundary,
                    boundary_provenance = d14.boundary_provenance,
                    post_d15_provenance = (;
                        source = d14.boundary_provenance.source,
                        stage = d14.boundary_provenance.stage,
                        route_status = d14.boundary_provenance.route_status,
                        current_peel_dimension = d14.boundary_provenance.current_peel_dimension,
                        last_completed_peel_dimension =
                            d14.boundary_provenance.last_completed_peel_dimension,
                        failure_code = d14.boundary_provenance.failure_code,
                        old_d15_boundary_cleared =
                            d14.boundary_provenance.old_d15_boundary_cleared,
                    ),
                    last_column_nonzero_count = d14.last_column_nonzero_count,
                    max_entry_term_count = d14.max_entry_term_count,
                ),
                :laurent_to_poly_case008_d14,
                :case008_d14_boundary_anchor,
            ),
        ],
    )
end

function laurent_to_poly_route_cases_by_id()
    return Dict(entry.id => entry for entry in laurent_to_poly_route_catalog().cases)
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
