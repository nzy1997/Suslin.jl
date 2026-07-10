using Test
using Oscar
using Suslin

include(joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl"))

function _route_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("route entry missing $(field)"))
    return getproperty(entry, field)
end

const _ROUTE_EXPECTATIONS = Dict(
    "laurent-to-poly-existing-normalization" => (;
        route = :existing_normalization,
        verifier_id = :laurent_to_poly_existing_normalization,
        source_case = :laurent_column_reduction_diagnostics,
        selected_entry_index = 1,
    ),
    "laurent-to-poly-general-ecp" => (;
        route = :general_ecp,
        verifier_id = :laurent_to_poly_general_ecp,
        source_case = :laurent_column_reduction_diagnostics,
        selected_entry_index = 2,
    ),
    "laurent-to-poly-case008-d14" => (;
        route = :case008_d14,
        verifier_id = :laurent_to_poly_case008_d14,
        source_case = "case_008",
        selected_entry_index = 1,
    ),
)

function _d14_fixture_from_route(entry)
    ring = entry.ring
    provenance = entry.provenance
    return (;
        case_id = provenance.source_case,
        source_case = provenance.source_case,
        source_cache_file = "case_008.jls",
        source_block = :column_transformation_upper_left_q_block,
        source_matrix_dimensions = (30, 30),
        source_column_transformation_dimensions = (60, 60),
        passed_peel_dimensions = (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15),
        first_failing_peel_dimension = 14,
        failing_column = entry.source_column,
        ring = ring.object,
        ring_description = ring.description,
        boundary_provenance = provenance.boundary_provenance,
        last_column_nonzero_count = 14,
        max_entry_term_count = 3734,
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

    diagnostic = if entry.route == :case008_d14
        ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(
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
    expected.stage === nothing || expected.stage in diagnostic.attempted_stages ||
        throw(ArgumentError("route $(entry.id) expected reducer stage is absent"))

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
    provenance.source_refs.laurent_noether.algorithm == :algorithm_6_3 ||
        throw(ArgumentError("route $(entry.id) has stale LaurentNoether provenance"))
    hasproperty(provenance, :route) || throw(ArgumentError("route $(entry.id) is missing route provenance"))
    provenance.route == entry.route || throw(ArgumentError("route $(entry.id) provenance route changed"))
    if entry.route == :case008_d14
        provenance.source_case == "case_008" ||
            throw(ArgumentError("route $(entry.id) case provenance changed"))
        provenance.source_fixture == :ToricBuilderCase008D14ColumnBoundary ||
            throw(ArgumentError("route $(entry.id) boundary fixture provenance changed"))
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
        merge(normalization, (; post_conversion_contract = merge(normalization.post_conversion_contract, (; selected_entry_must_be_polynomial_unit = false)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; post_conversion_contract = merge(normalization.post_conversion_contract, (; selected_entry_index = 2)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(d14, (; provenance = merge(d14.provenance, (; route = :stale_route)))),
    )
    @test_throws ArgumentError validate_laurent_to_poly_route_fixture(
        merge(normalization, (; source_column = [normalization.source_column[1] + one(normalization.ring.object); normalization.source_column[2:end]...])),
    )
end
