using Test
using Suslin
using Oscar

const LAURENT_FIXTURE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl")
const REQUIRED_LAURENT_FIXTURE_FIELDS = (
    :id,
    :kind,
    :ring_constructor,
    :ring,
    :dimensions,
    :inputs,
    :expected_relation,
    :provenance,
    :determinant_profile,
    :consumer_test_ids,
)

function _fixture_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _fixture_matrix_size(A)
    return (nrows(A), ncols(A))
end

function _fixture_determinant_classification(A)
    R = base_ring(A)
    d = det(A)
    d == one(R) && return "one"
    is_unit(d) && return "laurent_monomial_unit"
    return "non-unit"
end

function _assert_fixture_metadata(entry)
    for field in REQUIRED_LAURENT_FIXTURE_FIELDS
        _fixture_field(entry, field)
    end

    isempty(entry.consumer_test_ids) && throw(ArgumentError("fixture $(entry.id) must record at least one consumer test id"))
    hasproperty(entry.ring_constructor, :function_name) || throw(ArgumentError("fixture $(entry.id) missing ring constructor function name"))
    hasproperty(entry.ring_constructor, :coefficient) || throw(ArgumentError("fixture $(entry.id) missing ring constructor coefficient"))
    hasproperty(entry.ring_constructor, :variables) || throw(ArgumentError("fixture $(entry.id) missing ring constructor variables"))
    hasproperty(entry.ring, :description) || throw(ArgumentError("fixture $(entry.id) missing ring description"))
    hasproperty(entry.ring, :object) || throw(ArgumentError("fixture $(entry.id) missing ring object"))
    Suslin._require_laurent_polynomial_ring(entry.ring.object; label="fixture $(entry.id) ring")
    return true
end

function _assert_fixture_determinant_profile(entry)
    profile = entry.determinant_profile
    hasproperty(profile, :relevant) || throw(ArgumentError("fixture $(entry.id) missing determinant relevance"))
    profile.relevant || return true
    hasproperty(profile, :expected_class) || throw(ArgumentError("fixture $(entry.id) missing determinant class"))
    hasproperty(entry.inputs, :matrix) || throw(ArgumentError("fixture $(entry.id) determinant profile requires a matrix input"))

    actual_class = _fixture_determinant_classification(entry.inputs.matrix)
    actual_class == profile.expected_class || throw(ArgumentError("fixture $(entry.id) determinant classification $(profile.expected_class) does not match $(actual_class)"))
    return true
end

function _assert_solvable_linear_system_fixture(entry)
    A = entry.inputs.matrix
    rhs = entry.inputs.rhs
    solution = entry.inputs.expected_solution

    _fixture_matrix_size(A) == entry.dimensions.matrix || throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    _fixture_matrix_size(rhs) == entry.dimensions.rhs || throw(ArgumentError("fixture $(entry.id) rhs dimensions do not match metadata"))
    _fixture_matrix_size(solution) == entry.dimensions.solution || throw(ArgumentError("fixture $(entry.id) solution dimensions do not match metadata"))
    ncols(rhs) == 1 || throw(ArgumentError("fixture $(entry.id) rhs must be a column matrix"))
    ncols(solution) == 1 || throw(ArgumentError("fixture $(entry.id) solution must be a column matrix"))
    base_ring(A) == base_ring(rhs) || throw(ArgumentError("fixture $(entry.id) matrix and rhs use different rings"))
    base_ring(A) == base_ring(solution) || throw(ArgumentError("fixture $(entry.id) matrix and solution use different rings"))
    entry.expected_relation.kind == :linear_system_solution || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    A * solution == rhs || throw(ArgumentError("fixture $(entry.id) claimed solution does not satisfy matrix * solution == rhs"))
    return true
end

function _assert_unsolvable_linear_system_fixture(entry)
    A = entry.inputs.matrix
    rhs = entry.inputs.rhs
    certificate = entry.inputs.unsolvability_certificate

    _fixture_matrix_size(A) == entry.dimensions.matrix || throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    _fixture_matrix_size(rhs) == entry.dimensions.rhs || throw(ArgumentError("fixture $(entry.id) rhs dimensions do not match metadata"))
    base_ring(A) == base_ring(rhs) || throw(ArgumentError("fixture $(entry.id) matrix and rhs use different rings"))
    entry.expected_relation.kind == :linear_system_no_solution || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    certificate.kind == :zero_matrix_nonzero_rhs || throw(ArgumentError("fixture $(entry.id) uses an unsupported unsolvability certificate"))

    R = base_ring(A)
    A == zero_matrix(R, nrows(A), ncols(A)) || throw(ArgumentError("fixture $(entry.id) certificate requires a zero matrix"))
    i, j = certificate.rhs_index
    rhs[i, j] != zero(R) || throw(ArgumentError("fixture $(entry.id) certificate requires a nonzero rhs entry"))
    return true
end

function _assert_negative_normalization_fixture(entry)
    input_vector = entry.inputs.vector
    normalized_vector = entry.inputs.normalized_vector
    normalization_unit = entry.inputs.normalization_unit

    _fixture_matrix_size(input_vector) == entry.dimensions.vector || throw(ArgumentError("fixture $(entry.id) vector dimensions do not match metadata"))
    _fixture_matrix_size(normalized_vector) == entry.dimensions.normalized_vector || throw(ArgumentError("fixture $(entry.id) normalized dimensions do not match metadata"))
    base_ring(input_vector) == base_ring(normalized_vector) || throw(ArgumentError("fixture $(entry.id) vectors use different rings"))
    parent(normalization_unit) == base_ring(input_vector) || throw(ArgumentError("fixture $(entry.id) normalization unit uses the wrong ring"))
    entry.expected_relation.kind == :negative_exponent_normalization || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    normalization_unit * input_vector == normalized_vector || throw(ArgumentError("fixture $(entry.id) normalization relation failed"))
    return true
end

function _assert_toricbuilder_relation_fixture(entry)
    A = entry.inputs.matrix
    source = entry.inputs.source_matrix

    _fixture_matrix_size(A) == entry.dimensions.matrix || throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    _fixture_matrix_size(source) == entry.dimensions.source_matrix || throw(ArgumentError("fixture $(entry.id) source dimensions do not match metadata"))
    base_ring(A) == base_ring(source) || throw(ArgumentError("fixture $(entry.id) relation matrices use different rings"))
    entry.expected_relation.kind == :toricbuilder_left_inverse || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    entry.provenance.issue == "#19" || throw(ArgumentError("fixture $(entry.id) must trace to Issue 19"))
    entry.provenance.toricbuilder_commit == "fa7f82252d42fdc0b2726bc48af24ac4c70a8d73" || throw(ArgumentError("fixture $(entry.id) has unexpected ToricBuilder commit"))

    R = base_ring(A)
    source * A == identity_matrix(R, nrows(A)) || throw(ArgumentError("fixture $(entry.id) ToricBuilder relation failed"))
    return true
end

function validate_laurent_fixture(entry)
    _assert_fixture_metadata(entry)
    _assert_fixture_determinant_profile(entry)

    if entry.kind == :solvable_linear_system
        _assert_solvable_linear_system_fixture(entry)
    elseif entry.kind == :unsolvable_linear_system
        _assert_unsolvable_linear_system_fixture(entry)
    elseif entry.kind == :negative_exponent_normalization
        _assert_negative_normalization_fixture(entry)
    elseif entry.kind == :toricbuilder_relation
        _assert_toricbuilder_relation_fixture(entry)
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported kind $(entry.kind)"))
    end

    return true
end

function validate_laurent_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("Laurent fixture catalog missing cases"))
    isempty(catalog.cases) && throw(ArgumentError("Laurent fixture catalog must not be empty"))
    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) || throw(ArgumentError("Laurent fixture ids must be unique"))

    for entry in catalog.cases
        validate_laurent_fixture(entry)
    end
    return true
end

@testset "shared Laurent fixture catalog" begin
    @test isfile(LAURENT_FIXTURE_CATALOG_PATH)

    include(LAURENT_FIXTURE_CATALOG_PATH)
    catalog = LaurentFixtureCatalog.catalog()

    @test validate_laurent_fixture_catalog(catalog)

    fixture_by_id = Dict(entry.id => entry for entry in catalog.cases)
    @test haskey(fixture_by_id, "laurent-linear-system-solvable")
    @test haskey(fixture_by_id, "laurent-linear-system-unsolvable")
    @test haskey(fixture_by_id, "laurent-negative-exponent-normalization")
    @test haskey(fixture_by_id, "toricbuilder-factor-toric-block-3-pinv")

    solvable = fixture_by_id["laurent-linear-system-solvable"]
    bad_solution = copy(solvable.inputs.expected_solution)
    bad_solution[1, 1] += one(base_ring(bad_solution))
    bad_solvable = merge(solvable, (; inputs = merge(solvable.inputs, (; expected_solution = bad_solution))))
    @test_throws ArgumentError validate_laurent_fixture(bad_solvable)

    bad_det = merge(solvable, (; determinant_profile = merge(solvable.determinant_profile, (; expected_class = "non-unit"))))
    @test_throws ArgumentError validate_laurent_fixture(bad_det)

    toricbuilder = fixture_by_id["toricbuilder-factor-toric-block-3-pinv"]
    corrupted_source = copy(toricbuilder.inputs.source_matrix)
    corrupted_source[1, 1] += one(base_ring(corrupted_source))
    bad_toricbuilder = merge(toricbuilder, (; inputs = merge(toricbuilder.inputs, (; source_matrix = corrupted_source))))
    @test_throws ArgumentError validate_laurent_fixture(bad_toricbuilder)
end
