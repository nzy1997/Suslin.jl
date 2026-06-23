using Test
using Oscar
using Suslin

const PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

const REQUIRED_PARK_WOODBURN_POLYNOMIAL_IDS = Set([
    "pw-poly-univariate-sl3-fast-local-qq",
    "pw-poly-univariate-sln-disjoint-blocks-qq",
    "pw-poly-recursive-column-peel-gf2",
    "quillen-patched-substitution-witness-qq",
])

const REQUIRED_PARK_WOODBURN_POLYNOMIAL_NEGATIVE_IDS = Set([
    "pw-poly-det-not-one-control",
    "pw-poly-det-one-outside-witness-control",
    "pw-poly-wrong-route-control",
])

const SUPPORTED_PARK_WOODBURN_ROUTES = Set([
    :fast_local_sl3,
    :disjoint_local_blocks,
])

const STAGED_PARK_WOODBURN_ROUTES = Set([
    :recursive_column_peel,
    :quillen_patch,
    :quillen_patched_substitution,
])

function _pw_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _pw_require_matrix_over(matrix_value, R, n::Int, label)
    nrows(matrix_value) == n && ncols(matrix_value) == n ||
        throw(ArgumentError("$(label) must be a square matrix of fixture size"))
    base_ring(matrix_value) == R ||
        throw(ArgumentError("$(label) must be defined over the fixture ring"))
    return true
end

function _pw_assert_metadata(entry)
    for field in (
        :id,
        :role,
        :route,
        :status,
        :ring_constructor,
        :ring,
        :matrix,
        :determinant_expectation,
        :source_refs,
        :consumer_issue_ids,
    )
        _pw_field(entry, field)
    end

    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    entry.role isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) role must be a symbol"))

    route = _pw_field(entry, :route)
    status = _pw_field(entry, :status)
    route isa Symbol || throw(ArgumentError("fixture $(entry.id) route must be a symbol"))
    status isa Symbol || throw(ArgumentError("fixture $(entry.id) status must be a symbol"))

    if route in SUPPORTED_PARK_WOODBURN_ROUTES
        status == :supported ||
            throw(ArgumentError("fixture $(entry.id) with supported route $(route) must have status :supported"))
    elseif route in STAGED_PARK_WOODBURN_ROUTES
        status in (:supported, :staged, :blocked) &&
            status != :supported ||
            throw(ArgumentError("fixture $(entry.id) with route $(route) must not have status :supported"))
        if route in (:quillen_patch, :quillen_patched_substitution)
            consumer_issue_ids = _pw_field(entry, :consumer_issue_ids)
            "#105" in consumer_issue_ids ||
                throw(ArgumentError("quillen route fixture $(entry.id) must include issue #105"))
        end
    else
        throw(ArgumentError("fixture $(entry.id) uses unsupported route $(route)"))
    end

    ring_constructor = _pw_field(entry, :ring_constructor)
    _pw_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use polynomial_ring constructor metadata"))
    _pw_field(ring_constructor, :coefficient) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring constructor coefficient must be a string"))
    _pw_field(ring_constructor, :variables) isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be a tuple"))

    ring = _pw_field(entry, :ring)
    _pw_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    R = _pw_field(ring, :object)
    R isa MPolyRing ||
        throw(ArgumentError("fixture $(entry.id) must use an ordinary polynomial matrix ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("fixture $(entry.id) must use an exact coefficient ring"))
    generator_names = _pw_field(ring, :generator_names)
    generators = _pw_field(ring, :generators)
    generator_names isa Tuple && generators isa Tuple && length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring metadata is inconsistent"))
    ring_constructor.variables == generator_names ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables do not match ring metadata"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator parent mismatch"))

    matrix = _pw_field(entry, :matrix)
    size = nrows(matrix)
    matrix isa AbstractAlgebra.MatElem ||
        throw(ArgumentError("fixture $(entry.id) matrix must be a matrix"))
    size >= 3 || throw(ArgumentError("fixture $(entry.id) matrix size must be at least 3"))
    _pw_require_matrix_over(matrix, R, size, "matrix")
    all(idx -> matrix[idx, idx] isa Oscar.MPolyElem, 1:size) ||
        throw(ArgumentError("fixture $(entry.id) matrix entries must come from its ring"))

    determinant_expectation = _pw_field(entry, :determinant_expectation)
    determinant_expectation in (:one, :not_one) ||
        throw(ArgumentError("fixture $(entry.id) determinant expectation must be :one or :not_one"))
    if determinant_expectation == :one
        det(matrix) == one(R) ||
            throw(ArgumentError("fixture $(entry.id) matrix determinant must be one"))
    end

    entry.source_refs isa Tuple && !isempty(entry.source_refs) ||
        throw(ArgumentError("fixture $(entry.id) must include source refs"))
    entry.consumer_issue_ids isa Tuple && !isempty(entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) must include consumer issue ids"))
    all(id -> id isa AbstractString && startswith(id, "#"), entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must look like issue references"))
    return true
end

function validate_park_woodburn_polynomial_fixture(entry)
    _pw_assert_metadata(entry)
    return true
end

function validate_park_woodburn_polynomial_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain valid cases"))

    case_ids = [entry.id for entry in catalog.cases]
    length(case_ids) == length(unique(case_ids)) ||
        throw(ArgumentError("catalog valid case ids must be unique"))
    control_ids = [entry.id for entry in catalog.negative_controls]
    length(control_ids) == length(unique(control_ids)) ||
        throw(ArgumentError("catalog negative control ids must be unique"))

    for entry in catalog.cases
        validate_park_woodburn_polynomial_fixture(entry)
    end
    isempty(catalog.negative_controls) &&
        throw(ArgumentError("catalog must contain negative controls"))
    for entry in catalog.negative_controls
        try
            validate_park_woodburn_polynomial_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "Park-Woodburn polynomial fixture catalog" begin
    include(PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH)
    catalog = ParkWoodburnPolynomialFixtureCatalog.catalog()
    validate_park_woodburn_polynomial_fixture_catalog(catalog)
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()
    negatives = Dict(entry.id => entry for entry in catalog.negative_controls)

    @test REQUIRED_PARK_WOODBURN_POLYNOMIAL_IDS ⊆ Set(keys(entries))
    @test REQUIRED_PARK_WOODBURN_POLYNOMIAL_NEGATIVE_IDS ⊆ Set(keys(negatives))
    @test length(entries) >= 4

    for entry in values(negatives)
        @test_throws ArgumentError validate_park_woodburn_polynomial_fixture(entry)
    end

    staged_entry = entries["quillen-patched-substitution-witness-qq"]
    mutated_staged_route = merge(staged_entry, (; route = :fast_local_sl3))
    @test_throws ArgumentError validate_park_woodburn_polynomial_fixture(mutated_staged_route)

    supported_entry = entries["pw-poly-univariate-sl3-fast-local-qq"]
    mutated_supported_status = merge(supported_entry, (; status = :staged))
    @test_throws ArgumentError validate_park_woodburn_polynomial_fixture(mutated_supported_status)
end
