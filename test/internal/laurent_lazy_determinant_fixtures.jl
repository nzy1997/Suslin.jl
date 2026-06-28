using Test
using Suslin
using Oscar

const LAURENT_LAZY_DETERMINANT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "laurent_lazy_determinant_cases.jl")

const REQUIRED_LAZY_DETERMINANT_FIELDS = (
    :id,
    :kind,
    :ring,
    :dimensions,
    :inputs,
    :determinant_profile,
    :expected_correction,
    :normalizations,
    :negative_control,
    :provenance,
    :consumer_test_ids,
)

function _lazy_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("lazy determinant fixture $(get(entry, :id, "<unknown>")) missing field $(field)"))
    return getproperty(entry, field)
end

function _matrix_size(A)
    return (nrows(A), ncols(A))
end

function _fixture_supports_lazy_determinant_correction(entry)::Bool
    expected = _lazy_field(entry, :expected_correction)
    hasproperty(expected, :supported) || throw(ArgumentError("fixture $(entry.id) missing correction support flag"))
    hasproperty(entry, :inputs) && hasproperty(entry.inputs, :matrix) ||
        throw(ArgumentError("fixture $(entry.id) missing matrix input"))
    classification = classify_laurent_determinant(entry.inputs.matrix).classification
    actual_supported = classification in (:one, :laurent_monomial_unit, :permutation_sign_unit)
    return expected.supported && actual_supported
end

function _assert_negative_control_metadata(entry)
    negative_control = _lazy_field(entry, :negative_control)
    hasproperty(negative_control, :kind) || throw(ArgumentError("fixture $(entry.id) missing negative control kind"))
    hasproperty(negative_control, :base_case_id) || throw(ArgumentError("fixture $(entry.id) missing negative control base case id"))
    hasproperty(negative_control, :expected_failure) || throw(ArgumentError("fixture $(entry.id) missing negative control expected failure"))
    negative_control.base_case_id == entry.id ||
        throw(ArgumentError("fixture $(entry.id) negative control must refer to its base fixture id"))
    negative_control.kind == :metadata_mutation ||
        throw(ArgumentError("fixture $(entry.id) negative control kind changed"))
    negative_control.expected_failure == :determinant_class_metadata_mismatch ||
        throw(ArgumentError("fixture $(entry.id) negative control expected failure changed"))
    return true
end

function _assert_determinant_profile(entry)
    actual_profile = classify_laurent_determinant(entry.inputs.matrix)
    actual_profile.determinant == entry.determinant_profile.expected_determinant ||
        throw(ArgumentError("fixture $(entry.id) determinant does not match metadata"))
    actual_profile.classification == entry.determinant_profile.expected_class ||
        throw(ArgumentError("fixture $(entry.id) determinant class does not match metadata"))
    det(entry.inputs.matrix) == entry.determinant_profile.expected_determinant ||
        throw(ArgumentError("fixture $(entry.id) exact determinant does not match metadata"))
    entry.expected_correction.supported == (actual_profile.classification in (:one, :laurent_monomial_unit, :permutation_sign_unit)) ||
        throw(ArgumentError("fixture $(entry.id) supported-correction flag does not match determinant classification"))
    return actual_profile
end

function _assert_factorization_if_supported(core, label::AbstractString)
    nrows(core) >= 3 || return true
    factors = elementary_factorization(core)
    verify_factorization(core, factors) ||
        throw(ArgumentError("fixture $(label) elementary factorization did not verify"))
    return true
end

function _assert_supported_normalizations(entry)
    normalizations = _lazy_field(entry, :normalizations)
    hasproperty(normalizations, :row) || throw(ArgumentError("fixture $(entry.id) missing row normalization"))
    hasproperty(normalizations, :column) || throw(ArgumentError("fixture $(entry.id) missing column normalization"))

    row = normalizations.row
    column = normalizations.column
    hasproperty(row, :core) || throw(ArgumentError("fixture $(entry.id) missing row core"))
    hasproperty(row, :normalization) || throw(ArgumentError("fixture $(entry.id) missing row normalization object"))
    hasproperty(column, :core) || throw(ArgumentError("fixture $(entry.id) missing column core"))
    hasproperty(column, :correction_factor) || throw(ArgumentError("fixture $(entry.id) missing column correction factor"))

    verify_laurent_gl_normalization(entry.inputs.matrix, row.normalization) ||
        throw(ArgumentError("fixture $(entry.id) row normalization does not verify"))
    row.core == row.normalization.normalized_matrix ||
        throw(ArgumentError("fixture $(entry.id) row core must equal the normalized matrix"))
    entry.inputs.matrix * column.correction_factor == column.core ||
        throw(ArgumentError("fixture $(entry.id) column correction does not reconstruct the column core"))
    det(row.core) == one(base_ring(row.core)) ||
        throw(ArgumentError("fixture $(entry.id) row core must have determinant one"))
    det(column.core) == one(base_ring(column.core)) ||
        throw(ArgumentError("fixture $(entry.id) column core must have determinant one"))
    _assert_factorization_if_supported(row.core, "$(entry.id) row")
    _assert_factorization_if_supported(column.core, "$(entry.id) column")
    return true
end

function _assert_unsupported_normalizations(entry)
    entry.normalizations === nothing ||
        throw(ArgumentError("fixture $(entry.id) unsupported determinant must not provide normalizations"))
    !_fixture_supports_lazy_determinant_correction(entry) ||
        throw(ArgumentError("fixture $(entry.id) unsupported determinant must not report lazy correction support"))
    return true
end

function _assert_issue38_drift(entry)
    entry.id == "issue-38-q-block-lazy-determinant" || return true
    source = only(LaurentLazyDeterminantCases.ToricBuilderIssue38Cases.catalog().cases)
    u, v = source.ring.generators
    entry.inputs.matrix == source.inputs.matrix ||
        throw(ArgumentError("fixture $(entry.id) drifted from the Issue #38 source matrix"))
    entry.determinant_profile.expected_determinant == source.determinant_profile.expected_determinant ||
        throw(ArgumentError("fixture $(entry.id) drifted from the Issue #38 source determinant"))
    entry.determinant_profile.expected_determinant == u * v ||
        throw(ArgumentError("fixture $(entry.id) expected determinant must be u*v"))
    entry.normalizations.row.core == source.normalizations.row.core ||
        throw(ArgumentError("fixture $(entry.id) row core drifted from the Issue #38 source fixture"))
    entry.normalizations.column.core == source.normalizations.column.core ||
        throw(ArgumentError("fixture $(entry.id) column core drifted from the Issue #38 source fixture"))
    return true
end

function validate_laurent_lazy_determinant_fixture(entry)
    for field in REQUIRED_LAZY_DETERMINANT_FIELDS
        _lazy_field(entry, field)
    end
    isempty(entry.consumer_test_ids) && throw(ArgumentError("fixture $(entry.id) must record at least one consumer test id"))
    _matrix_size(entry.inputs.matrix) == entry.dimensions.matrix ||
        throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    _assert_negative_control_metadata(entry)
    actual_profile = _assert_determinant_profile(entry)
    if actual_profile.classification in (:one, :laurent_monomial_unit, :permutation_sign_unit)
        _assert_supported_normalizations(entry)
    else
        _assert_unsupported_normalizations(entry)
    end
    _assert_issue38_drift(entry)
    return true
end

function validate_laurent_lazy_determinant_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("lazy determinant catalog missing cases"))
    isempty(catalog.cases) && throw(ArgumentError("lazy determinant catalog must not be empty"))
    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) || throw(ArgumentError("lazy determinant fixture ids must be unique"))
    for entry in catalog.cases
        validate_laurent_lazy_determinant_fixture(entry)
    end
    return true
end

@testset "lazy Laurent determinant fixture catalog" begin
    @test isfile(LAURENT_LAZY_DETERMINANT_FIXTURE_PATH)

    include(LAURENT_LAZY_DETERMINANT_FIXTURE_PATH)
    catalog = LaurentLazyDeterminantCases.catalog()
    @test validate_laurent_lazy_determinant_catalog(catalog)

    fixture_by_id = Dict(entry.id => entry for entry in catalog.cases)
    @test haskey(fixture_by_id, "issue-38-q-block-lazy-determinant")
    @test haskey(fixture_by_id, "determinant-one-triangular")
    @test haskey(fixture_by_id, "monomial-unit-row-column-cores")
    @test haskey(fixture_by_id, "non-unit-determinant-negative")

    monomial = fixture_by_id["monomial-unit-row-column-cores"]
    bad_monomial = merge(
        monomial,
        (; determinant_profile = merge(monomial.determinant_profile, (; expected_class = :one))),
    )
    @test_throws ArgumentError validate_laurent_lazy_determinant_fixture(bad_monomial)

    non_unit = fixture_by_id["non-unit-determinant-negative"]
    @test !non_unit.expected_correction.supported
    @test !_fixture_supports_lazy_determinant_correction(non_unit)
end
