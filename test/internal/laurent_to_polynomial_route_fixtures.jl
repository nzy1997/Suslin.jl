using Test
using Oscar
using Suslin

include(joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl"))

const _D14Boundary = LaurentFixtureCatalog.ToricBuilderCase008D14ColumnBoundary

function _record_field(record, field::Symbol, label)
    hasproperty(record, field) || throw(ArgumentError("$(label) missing $(field)"))
    return getproperty(record, field)
end

function _route_field(entry, field::Symbol)
    return _record_field(entry, field, "route entry")
end

function _entry_term_count(entry)
    iszero(entry) && return 0
    return length(collect(coefficients(entry)))
end

function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) ||
        throw(ArgumentError("route diagnostic missing stage details"))
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _validate_stage_detail_shape(entry, diagnostic)
    hasproperty(diagnostic, :stage_details) ||
        throw(ArgumentError("route $(entry.id) reducer diagnostic has no stage details"))
    diagnostic.stage_details isa Tuple ||
        throw(ArgumentError("route $(entry.id) reducer stage details are not a tuple"))
    length(diagnostic.stage_details) == length(diagnostic.attempted_stages) ||
        throw(ArgumentError("route $(entry.id) reducer stage details are stale"))
    all(detail -> detail isa NamedTuple, diagnostic.stage_details) ||
        throw(ArgumentError("route $(entry.id) reducer stage details are malformed"))
end

function _validate_expected_stage_detail(entry, diagnostic, expected)
    expected.stage === nothing && return nothing
    detail = _diagnostic_stage_detail(diagnostic, expected.stage)
    detail !== nothing ||
        throw(ArgumentError("route $(entry.id) expected reducer stage is absent"))
    detail.outcome == expected.stage_outcome ||
        throw(ArgumentError("route $(entry.id) reducer stage outcome changed"))

    if expected.stage_outcome == :delegated_to_polynomial
        detail.normalized_status == expected.normalized_status ||
            throw(ArgumentError("route $(entry.id) normalized reducer status changed"))
        detail.normalized_failure_code == expected.normalized_failure_code ||
            throw(ArgumentError("route $(entry.id) normalized reducer failure code changed"))
    elseif expected.stage_outcome == :staged_boundary
        detail.boundary == expected.boundary ||
            throw(ArgumentError("route $(entry.id) reducer boundary changed"))
        detail.requires_descent_measure == expected.requires_descent_measure ||
            throw(ArgumentError("route $(entry.id) descent-measure requirement changed"))
        detail.certified_descent_scope == expected.certified_descent_scope ||
            throw(ArgumentError("route $(entry.id) certified descent scope changed"))
        detail.next_boundary == expected.next_boundary ||
            throw(ArgumentError("route $(entry.id) next boundary changed"))
        detail.requires_link_witness == expected.requires_link_witness ||
            throw(ArgumentError("route $(entry.id) link-witness requirement changed"))
        detail.requires_endpoint_reduction == expected.requires_endpoint_reduction ||
            throw(ArgumentError("route $(entry.id) endpoint-reduction requirement changed"))
        detail.requires_laurent_normality_replay == expected.requires_laurent_normality_replay ||
            throw(ArgumentError("route $(entry.id) Laurent normality replay requirement changed"))
        detail.requires_recursive_peel_integration ==
            expected.requires_recursive_peel_integration ||
            throw(ArgumentError("route $(entry.id) recursive peel requirement changed"))
        detail.fallback_policy == expected.fallback_policy ||
            throw(ArgumentError("route $(entry.id) boundary fallback policy changed"))
    else
        throw(ArgumentError("route $(entry.id) has unsupported expected stage outcome"))
    end
    return detail
end

const _ROUTE_EXPECTATIONS = Dict(
    "laurent-to-poly-existing-normalization" => (;
        route = :existing_normalization,
        verifier_id = :laurent_to_poly_existing_normalization,
        source_case = :laurent_column_reduction_diagnostics,
        selected_entry_index = 1,
        selected_entry_role = :existing_normalization_anchor,
    ),
    "laurent-to-poly-general-ecp" => (;
        route = :general_ecp,
        verifier_id = :laurent_to_poly_general_ecp,
        source_case = :laurent_column_reduction_diagnostics,
        selected_entry_index = 2,
        selected_entry_role = :general_ecp_anchor,
    ),
    "laurent-to-poly-case008-d14" => (;
        route = :case008_d14,
        verifier_id = :laurent_to_poly_case008_d14,
        source_case = "case_008",
        selected_entry_index = 1,
        selected_entry_role = :case008_d14_boundary_anchor,
    ),
)

function _d14_fixture_from_route(entry)
    ring = entry.ring
    provenance = entry.provenance
    label = "route $(entry.id) d14 provenance"
    _record_field(provenance, :dimension, label) ==
        _record_field(provenance, :first_failing_peel_dimension, label) ||
        throw(ArgumentError("route $(entry.id) d14 dimension metadata is inconsistent"))
    return (;
        case_id = _record_field(provenance, :case_id, label),
        source_case = _record_field(provenance, :source_case, label),
        source_cache_file = _record_field(provenance, :source_cache_file, label),
        source_block = _record_field(provenance, :source_block, label),
        source_matrix_dimensions = _record_field(provenance, :source_matrix_dimensions, label),
        source_column_transformation_dimensions =
            _record_field(provenance, :source_column_transformation_dimensions, label),
        passed_peel_dimensions = _record_field(provenance, :passed_peel_dimensions, label),
        first_failing_peel_dimension = _record_field(provenance, :dimension, label),
        failing_column = entry.source_column,
        ring = ring.object,
        ring_description = ring.description,
        boundary_provenance = _record_field(provenance, :boundary_provenance, label),
        last_column_nonzero_count = _record_field(provenance, :last_column_nonzero_count, label),
        max_entry_term_count = _record_field(provenance, :max_entry_term_count, label),
    )
end

function validate_laurent_to_poly_route_fixture(entry)
    ring = _route_field(entry, :ring)
    column = _route_field(entry, :source_column)
    constructor = _route_field(entry, :ring_constructor)
    provenance = _route_field(entry, :provenance)
    expected = _route_field(entry, :expected_reducer)
    contract = _route_field(entry, :post_conversion_contract)
    expectation = get(_ROUTE_EXPECTATIONS, entry.id, nothing)
    expectation === nothing && throw(ArgumentError("unknown route $(entry.id)"))

    parent(first(column)) == ring.object ||
        throw(ArgumentError("route $(entry.id) column has the wrong ring"))
    length(ring.generators) == 2 ||
        throw(ArgumentError("route $(entry.id) does not have two Laurent generators"))
    base_ring(ring.object) == GF(2) ||
        throw(ArgumentError("route $(entry.id) does not have GF(2) coefficients"))
    ring.description == "GF(2)[$(ring.generators[1])^+/-1, $(ring.generators[2])^+/-1]" ||
        throw(ArgumentError("route $(entry.id) ring description is stale"))
    constructor.function_name == :suslin_laurent_polynomial_ring ||
        throw(ArgumentError("route $(entry.id) has the wrong ring constructor"))
    Tuple(string.(gens(ring.object))) == Tuple(string.(ring.generators)) ||
        throw(ArgumentError("route $(entry.id) Laurent generators changed"))
    constructor.coefficient == "GF(2)" ||
        throw(ArgumentError("route $(entry.id) constructor coefficient changed"))
    constructor.variables == Tuple(string.(gens(ring.object))) ||
        throw(ArgumentError("route $(entry.id) constructor variables do not match the ring"))
    constructor.variables == Tuple(string.(ring.generators)) ||
        throw(ArgumentError("route $(entry.id) constructor variables do not match the route ring"))

    Suslin._require_laurent_polynomial_ring(ring.object; label="route $(entry.id) ring")
    Suslin._laurent_descent_column_support_fingerprint(column) == entry.source_fingerprint ||
        throw(ArgumentError("route $(entry.id) source fingerprint is stale"))
    1 <= entry.selected_entry_index <= length(column) ||
        throw(ArgumentError("route $(entry.id) selected entry is out of bounds"))
    entry.selected_entry_index == expectation.selected_entry_index ||
        throw(ArgumentError("route $(entry.id) selected entry changed"))
    selected_entry = column[entry.selected_entry_index]

    diagnostic = if entry.route == :case008_d14
        _D14Boundary.validate_boundary_fixture(
            _d14_fixture_from_route(entry),
        ) == :ok || throw(ArgumentError("route $(entry.id) boundary fixture is stale or non-unimodular"))
        Suslin.diagnose_unimodular_column_reduction(
            column,
            ring.object;
            assume_unimodular = true,
            laurent_large_support_diagnostic_decline = true,
        )
    else
        Suslin.is_unimodular_column(column, ring.object) ||
            throw(ArgumentError("route $(entry.id) source column is not unimodular"))
        Suslin.diagnose_unimodular_column_reduction(column, ring.object)
    end
    diagnostic.status == expected.status ||
        throw(ArgumentError("route $(entry.id) reducer status changed"))
    diagnostic.failure_code == expected.failure_code ||
        throw(ArgumentError("route $(entry.id) reducer failure code changed"))
    _validate_stage_detail_shape(entry, diagnostic)
    _validate_expected_stage_detail(entry, diagnostic, expected)

    entry.route == expectation.route || throw(ArgumentError("route $(entry.id) identifier changed"))
    entry.verifier_id == expectation.verifier_id ||
        throw(ArgumentError("route $(entry.id) verifier identifier changed"))
    hasproperty(provenance, :source_refs) ||
        throw(ArgumentError("route $(entry.id) is missing source provenance"))
    provenance.source_case == expectation.source_case ||
        throw(ArgumentError("route $(entry.id) source provenance changed"))
    provenance.source_refs.laurent_to_poly.author == "Park" ||
        throw(ArgumentError("route $(entry.id) LaurentToPoly author provenance changed"))
    provenance.source_refs.laurent_to_poly.algorithm == :algorithm_6_1 ||
        throw(ArgumentError("route $(entry.id) has stale LaurentToPoly provenance"))
    provenance.source_refs.laurent_to_poly.name == "LaurentToPoly" ||
        throw(ArgumentError("route $(entry.id) LaurentToPoly name provenance changed"))
    provenance.source_refs.laurent_to_poly.role == :laurent_to_polynomial_conversion ||
        throw(ArgumentError("route $(entry.id) LaurentToPoly role provenance changed"))
    provenance.source_refs.laurent_noether.author == "Park" ||
        throw(ArgumentError("route $(entry.id) LaurentNoether author provenance changed"))
    provenance.source_refs.laurent_noether.algorithm == :algorithm_6_3 ||
        throw(ArgumentError("route $(entry.id) has stale LaurentNoether provenance"))
    provenance.source_refs.laurent_noether.name == "LaurentNoether" ||
        throw(ArgumentError("route $(entry.id) LaurentNoether name provenance changed"))
    provenance.source_refs.laurent_noether.role == :laurent_variable_change_normalization ||
        throw(ArgumentError("route $(entry.id) LaurentNoether role provenance changed"))
    hasproperty(provenance, :route) || throw(ArgumentError("route $(entry.id) is missing route provenance"))
    provenance.route == entry.route || throw(ArgumentError("route $(entry.id) provenance route changed"))
    if entry.route == :case008_d14
        canonical_d14 = _D14Boundary.boundary_fixture()
        provenance.source_case == "case_008" ||
            throw(ArgumentError("route $(entry.id) case provenance changed"))
        provenance.source_fixture == :ToricBuilderCase008D14ColumnBoundary ||
            throw(ArgumentError("route $(entry.id) boundary fixture provenance changed"))
        provenance.case_id == canonical_d14.case_id ||
            throw(ArgumentError("route $(entry.id) d14 case id changed"))
        provenance.dimension == canonical_d14.first_failing_peel_dimension ||
            throw(ArgumentError("route $(entry.id) d14 dimension changed"))
        provenance.dimension == 14 ||
            throw(ArgumentError("route $(entry.id) d14 dimension must be 14"))
        provenance.source_cache_file == canonical_d14.source_cache_file ||
            throw(ArgumentError("route $(entry.id) d14 cache provenance changed"))
        provenance.source_block == canonical_d14.source_block ||
            throw(ArgumentError("route $(entry.id) d14 source block changed"))
        provenance.source_matrix_dimensions == canonical_d14.source_matrix_dimensions ||
            throw(ArgumentError("route $(entry.id) d14 source matrix dimensions changed"))
        provenance.source_column_transformation_dimensions ==
            canonical_d14.source_column_transformation_dimensions ||
            throw(ArgumentError("route $(entry.id) d14 transformation dimensions changed"))
        provenance.passed_peel_dimensions == canonical_d14.passed_peel_dimensions ||
            throw(ArgumentError("route $(entry.id) d14 post-d15 peel provenance changed"))
        provenance.first_failing_peel_dimension == canonical_d14.first_failing_peel_dimension ||
            throw(ArgumentError("route $(entry.id) d14 first failing dimension changed"))
        provenance.boundary_status == :current_staged_d14_boundary ||
            throw(ArgumentError("route $(entry.id) d14 boundary status changed"))
        provenance.boundary_provenance == canonical_d14.boundary_provenance ||
            throw(ArgumentError("route $(entry.id) d14 boundary provenance changed"))
        provenance.post_d15_provenance == (;
            source = canonical_d14.boundary_provenance.source,
            stage = canonical_d14.boundary_provenance.stage,
            route_status = canonical_d14.boundary_provenance.route_status,
            current_peel_dimension = canonical_d14.boundary_provenance.current_peel_dimension,
            last_completed_peel_dimension =
                canonical_d14.boundary_provenance.last_completed_peel_dimension,
            failure_code = canonical_d14.boundary_provenance.failure_code,
            old_d15_boundary_cleared =
                canonical_d14.boundary_provenance.old_d15_boundary_cleared,
        ) || throw(ArgumentError("route $(entry.id) d14 post-d15 provenance changed"))
    end
    contract.preserves_unimodularity === true ||
        throw(ArgumentError("route $(entry.id) does not preserve unimodularity"))
    contract.polynomial_target === true ||
        throw(ArgumentError("route $(entry.id) has no polynomial target contract"))
    contract.selected_entry_must_be_polynomial_unit === true ||
        throw(ArgumentError("route $(entry.id) selected entry contract changed"))
    hasproperty(contract, :selected_entry_index) ||
        throw(ArgumentError("route $(entry.id) has no selected entry contract"))
    contract.selected_entry_index == entry.selected_entry_index ||
        throw(ArgumentError("route $(entry.id) selected entry contract is stale"))
    contract.selected_entry_role == expectation.selected_entry_role ||
        throw(ArgumentError("route $(entry.id) selected entry role changed"))
    contract.selected_source_fingerprint ==
        Suslin._laurent_descent_column_support_fingerprint([selected_entry]) ||
        throw(ArgumentError("route $(entry.id) selected entry fingerprint is stale"))
    contract.selected_source_term_count == _entry_term_count(selected_entry) ||
        throw(ArgumentError("route $(entry.id) selected entry term count is stale"))
    contract.selected_source_is_unit === is_unit(selected_entry) ||
        throw(ArgumentError("route $(entry.id) selected entry unit profile is stale"))
    entry.consumer_test_ids == (
        "issue-351-laurent-to-poly-route-fixtures",
        "issue-351-laurent-to-poly-route-consumer",
    ) || throw(ArgumentError("route $(entry.id) downstream consumers changed"))
    return :ok
end

@testset "LaurentToPoly route fixture catalog" begin
    catalog = LaurentFixtureCatalog.laurent_to_poly_route_catalog()
    @test Set(entry.id for entry in catalog.cases) == Set([
        "laurent-to-poly-existing-normalization",
        "laurent-to-poly-general-ecp",
        "laurent-to-poly-case008-d14",
    ])
    @test length(catalog.cases) == 3
    by_id = LaurentFixtureCatalog.laurent_to_poly_route_cases_by_id()
    @test Set(keys(by_id)) == Set(entry.id for entry in catalog.cases)
    for entry in catalog.cases
        @test validate_laurent_to_poly_route_fixture(entry) == :ok
    end

    d14 = by_id["laurent-to-poly-case008-d14"]
    d14_nonunimodular_column = [
        (gens(d14.ring.object)[2] + one(d14.ring.object)) * value
        for value in d14.source_column
    ]
    d14_nonunimodular = merge(d14, (;
        source_column = d14_nonunimodular_column,
        source_fingerprint = Suslin._laurent_descent_column_support_fingerprint(d14_nonunimodular_column),
    ))
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(d14_nonunimodular)

    normalization = by_id["laurent-to-poly-existing-normalization"]
    x, y = gens(normalization.ring.object)
    swapped = merge(normalization, (; ring = merge(normalization.ring, (; generators = (y, x)))))
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(swapped)
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; ring_constructor = merge(normalization.ring_constructor, (; coefficient = "GF(3)")))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; ring_constructor = merge(normalization.ring_constructor, (; variables = ("y", "x"))))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; selected_entry_index = length(normalization.source_column) + 1)),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; selected_entry_index = 2)),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(
            normalization,
            (;
                expected_reducer = merge(
                    normalization.expected_reducer,
                    (; stage_outcome = :supported),
                ),
            ),
        ),
    )
    general = by_id["laurent-to-poly-general-ecp"]
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(
            general,
            (;
                expected_reducer = merge(
                    general.expected_reducer,
                    (; requires_descent_measure = false),
                ),
            ),
        ),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; post_conversion_contract = merge(normalization.post_conversion_contract, (; selected_entry_must_be_polynomial_unit = false)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; post_conversion_contract = merge(normalization.post_conversion_contract, (; selected_entry_index = 2)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; post_conversion_contract = merge(normalization.post_conversion_contract, (; selected_entry_role = :wrong_anchor)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(
            normalization,
            (;
                post_conversion_contract = merge(
                    normalization.post_conversion_contract,
                    (;
                        selected_source_fingerprint =
                            Suslin._laurent_descent_column_support_fingerprint([normalization.source_column[2]]),
                    ),
                ),
            ),
        ),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; post_conversion_contract = merge(normalization.post_conversion_contract, (; selected_source_term_count = 0)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(
            normalization,
            (;
                post_conversion_contract = merge(
                    normalization.post_conversion_contract,
                    (; selected_source_is_unit = !normalization.post_conversion_contract.selected_source_is_unit),
                ),
            ),
        ),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(d14, (; provenance = merge(d14.provenance, (; route = :stale_route)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(d14, (; provenance = merge(d14.provenance, (; dimension = 15)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(
            d14,
            (;
                provenance = merge(
                    d14.provenance,
                    (;
                        post_d15_provenance = merge(
                            d14.provenance.post_d15_provenance,
                            (; old_d15_boundary_cleared = false),
                        ),
                    ),
                ),
            ),
        ),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(
            normalization,
            (;
                provenance = merge(
                    normalization.provenance,
                    (;
                        source_refs = merge(
                            normalization.provenance.source_refs,
                            (;
                                laurent_noether = merge(
                                    normalization.provenance.source_refs.laurent_noether,
                                    (; author = "stale"),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; source_column = [normalization.source_column[1] + one(normalization.ring.object); normalization.source_column[2:end]...])),
    )
end
