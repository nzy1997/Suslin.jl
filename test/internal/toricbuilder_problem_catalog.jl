using Test
using Suslin
using Oscar

const TORICBUILDER_PROBLEM_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_laurent_problem_catalog.jl")

const REQUIRED_TORICBUILDER_PROBLEM_FIELDS = (
    :id,
    :kind,
    :source_fixture,
    :ring,
    :dimensions,
    :matrix,
    :determinant_profile,
    :expected_current_status,
    :verifier,
    :provenance,
    :consumers,
)

const REQUIRED_TORICBUILDER_PROBLEM_IDS = Set([
    "toricbuilder-issue-38-q-block",
    "toricbuilder-factor-toric-block-3-qinv",
    "toricbuilder-factor-toric-block-3-pinv",
    "laurent-block-local-40x40",
    "laurent-block-local-48x48",
])

const ALLOWED_TORICBUILDER_PROBLEM_STATUSES = Set([
    :unsupported_now,
    :supported_column_peel,
    :verified_contract,
    :supported_block_local,
    :target_acceptance,
])

function _problem_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("ToricBuilder problem entry missing field $(field)"))
    return getproperty(entry, field)
end

function _problem_matrix_size(A)
    return (nrows(A), ncols(A))
end

function _problem_determinant_classification(A)
    profile = Suslin.classify_laurent_determinant(A)
    return profile.classification
end

function _assert_nonempty_metadata(value, label::AbstractString)
    value === nothing && throw(ArgumentError("$(label) must not be nothing"))
    if value isa AbstractString
        isempty(value) && throw(ArgumentError("$(label) must not be empty"))
    elseif value isa Tuple || value isa AbstractVector
        isempty(value) && throw(ArgumentError("$(label) must not be empty"))
    elseif value isa NamedTuple
        isempty(keys(value)) && throw(ArgumentError("$(label) must not be empty"))
    end
    return true
end

function validate_toricbuilder_problem_entry(entry)
    for field in REQUIRED_TORICBUILDER_PROBLEM_FIELDS
        _problem_field(entry, field)
    end

    _assert_nonempty_metadata(entry.id, "problem id")
    _assert_nonempty_metadata(entry.provenance, "problem $(entry.id) provenance")
    _assert_nonempty_metadata(entry.consumers, "problem $(entry.id) consumers")
    hasproperty(entry.consumers, :milestone) ||
        throw(ArgumentError("problem $(entry.id) missing consumer milestone"))
    hasproperty(entry.consumers, :issues) ||
        throw(ArgumentError("problem $(entry.id) missing consumer issues"))
    hasproperty(entry.consumers, :tests) ||
        throw(ArgumentError("problem $(entry.id) missing consumer tests"))
    isempty(entry.consumers.issues) &&
        throw(ArgumentError("problem $(entry.id) must record at least one consumer issue"))
    isempty(entry.consumers.tests) &&
        throw(ArgumentError("problem $(entry.id) must record at least one consumer test"))
    entry.expected_current_status in ALLOWED_TORICBUILDER_PROBLEM_STATUSES ||
        throw(ArgumentError("problem $(entry.id) has unsupported expected status $(entry.expected_current_status)"))
    if entry.expected_current_status == :supported_column_peel
        hasproperty(entry, :expected_suslin_path) ||
            throw(ArgumentError("problem $(entry.id) missing expected Suslin path"))
        entry.expected_suslin_path == :laurent_column_peel ||
            throw(ArgumentError("problem $(entry.id) expected Laurent column-peel path"))
    end

    hasproperty(entry.verifier, :path) ||
        throw(ArgumentError("problem $(entry.id) missing verifier path"))
    isfile(joinpath(@__DIR__, "..", entry.verifier.path)) ||
        throw(ArgumentError("problem $(entry.id) verifier path $(entry.verifier.path) does not exist"))

    hasproperty(entry.dimensions, :matrix) ||
        throw(ArgumentError("problem $(entry.id) missing matrix dimensions"))
    _problem_matrix_size(entry.matrix) == entry.dimensions.matrix ||
        throw(ArgumentError("problem $(entry.id) matrix dimensions do not match metadata"))

    hasproperty(entry.determinant_profile, :expected_class) ||
        throw(ArgumentError("problem $(entry.id) missing determinant class"))
    hasproperty(entry.determinant_profile, :expected_determinant) ||
        throw(ArgumentError("problem $(entry.id) missing expected determinant"))
    entry.determinant_profile.expected_determinant === nothing &&
        throw(ArgumentError("problem $(entry.id) expected determinant must not be nothing"))
    _problem_determinant_classification(entry.matrix) == entry.determinant_profile.expected_class ||
        throw(ArgumentError("problem $(entry.id) determinant class does not match metadata"))
    det(entry.matrix) == entry.determinant_profile.expected_determinant ||
        throw(ArgumentError("problem $(entry.id) determinant does not match metadata"))

    return true
end

function validate_toricbuilder_problem_catalog(catalog)
    hasproperty(catalog, :cases) ||
        throw(ArgumentError("ToricBuilder problem catalog missing cases"))
    isempty(catalog.cases) &&
        throw(ArgumentError("ToricBuilder problem catalog must not be empty"))

    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) ||
        throw(ArgumentError("ToricBuilder problem catalog ids must be unique"))
    REQUIRED_TORICBUILDER_PROBLEM_IDS ⊆ Set(ids) ||
        throw(ArgumentError("ToricBuilder problem catalog missing required IDs"))

    for entry in catalog.cases
        validate_toricbuilder_problem_entry(entry)
    end

    return true
end

@testset "ToricBuilder Laurent problem catalog" begin
    @test isfile(TORICBUILDER_PROBLEM_CATALOG_PATH)

    include(TORICBUILDER_PROBLEM_CATALOG_PATH)
    catalog = ToricBuilderLaurentProblemCatalog.catalog()
    @test validate_toricbuilder_problem_catalog(catalog)

    by_id = Dict(entry.id => entry for entry in catalog.cases)
    @test haskey(by_id, "toricbuilder-issue-38-q-block")
    @test haskey(by_id, "toricbuilder-factor-toric-block-3-qinv")
    @test haskey(by_id, "toricbuilder-factor-toric-block-3-pinv")
    @test haskey(by_id, "laurent-block-local-40x40")
    @test haskey(by_id, "laurent-block-local-48x48")
    @test by_id["toricbuilder-issue-38-q-block"].expected_current_status == :unsupported_now
    @test by_id["toricbuilder-factor-toric-block-3-qinv"].expected_current_status == :supported_column_peel
    @test by_id["toricbuilder-factor-toric-block-3-qinv"].expected_suslin_path == :laurent_column_peel
    @test by_id["toricbuilder-factor-toric-block-3-pinv"].expected_current_status == :supported_column_peel
    @test by_id["toricbuilder-factor-toric-block-3-pinv"].expected_suslin_path == :laurent_column_peel
    @test by_id["laurent-block-local-40x40"].expected_current_status == :supported_block_local
    @test by_id["laurent-block-local-48x48"].expected_current_status == :target_acceptance

    duplicate_cases = [catalog.cases; (merge(first(catalog.cases), (; id = catalog.cases[2].id,)))]
    @test_throws ArgumentError validate_toricbuilder_problem_catalog((; cases = duplicate_cases))

    missing_provenance = merge(first(catalog.cases), (; provenance = (;),))
    @test_throws ArgumentError validate_toricbuilder_problem_entry(missing_provenance)

    bad_determinant = merge(
        first(catalog.cases),
        (;
            determinant_profile = merge(
                first(catalog.cases).determinant_profile,
                (; expected_determinant = one(base_ring(first(catalog.cases).matrix))),
            ),
        ),
    )
    @test_throws ArgumentError validate_toricbuilder_problem_entry(bad_determinant)
end
