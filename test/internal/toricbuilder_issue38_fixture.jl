using Test
using Suslin
using Oscar

const TORICBUILDER_ISSUE38_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_issue38_cases.jl")

const REQUIRED_TORICBUILDER_ISSUE38_FIELDS = (
    :id,
    :kind,
    :ring,
    :dimensions,
    :inputs,
    :determinant_profile,
    :normalizations,
    :expected_current_status,
    :provenance,
    :consumer_test_ids,
)

function _fixture_matrix_size(A)
    return (nrows(A), ncols(A))
end

function _require_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _assert_issue38_metadata(entry)
    for field in REQUIRED_TORICBUILDER_ISSUE38_FIELDS
        _require_field(entry, field)
    end

    isempty(entry.consumer_test_ids) && throw(ArgumentError("fixture $(entry.id) must record at least one consumer test id"))
    entry.kind == :toricbuilder_q_block || throw(ArgumentError("fixture $(entry.id) has unexpected kind $(entry.kind)"))
    hasproperty(entry.ring, :description) || throw(ArgumentError("fixture $(entry.id) missing ring description"))
    hasproperty(entry.ring, :object) || throw(ArgumentError("fixture $(entry.id) missing ring object"))
    hasproperty(entry.ring, :generators) || throw(ArgumentError("fixture $(entry.id) missing ring generators"))
    hasproperty(entry.ring, :variables) || throw(ArgumentError("fixture $(entry.id) missing ring variables"))
    hasproperty(entry.provenance, :issue) || throw(ArgumentError("fixture $(entry.id) missing provenance issue"))
    hasproperty(entry.provenance, :reported_main_commit) || throw(ArgumentError("fixture $(entry.id) missing provenance commit"))
    hasproperty(entry.inputs, :matrix) || throw(ArgumentError("fixture $(entry.id) missing matrix input"))
    hasproperty(entry.dimensions, :matrix) || throw(ArgumentError("fixture $(entry.id) missing matrix dimensions"))

    _fixture_matrix_size(entry.inputs.matrix) == (6, 6) || throw(ArgumentError("fixture $(entry.id) must be a 6x6 matrix"))
    entry.dimensions.matrix == (6, 6) || throw(ArgumentError("fixture $(entry.id) matrix dimensions must be 6x6"))
    return true
end

function _assert_issue38_determinant_profile(entry)
    profile = _require_field(entry, :determinant_profile)
    hasproperty(profile, :expected_determinant) || throw(ArgumentError("fixture $(entry.id) missing expected determinant"))
    hasproperty(profile, :expected_class) || throw(ArgumentError("fixture $(entry.id) missing expected determinant class"))

    A = entry.inputs.matrix
    actual_profile = classify_laurent_determinant(A)
    det(A) == profile.expected_determinant || throw(ArgumentError("fixture $(entry.id) determinant does not match metadata"))
    actual_profile.classification == profile.expected_class || throw(ArgumentError("fixture $(entry.id) determinant class does not match metadata"))
    return true
end

function _assert_issue38_row_normalization(entry)
    normalizations = _require_field(entry, :normalizations)
    hasproperty(normalizations, :row) || throw(ArgumentError("fixture $(entry.id) missing row normalization metadata"))

    row = normalizations.row
    hasproperty(row, :core) || throw(ArgumentError("fixture $(entry.id) missing row normalization core"))
    hasproperty(row, :normalization) || throw(ArgumentError("fixture $(entry.id) missing row normalization object"))
    det(row.core) == one(base_ring(row.core)) || throw(ArgumentError("fixture $(entry.id) row normalization core must have determinant one"))
    verify_laurent_gl_normalization(entry.inputs.matrix, row.normalization) || throw(ArgumentError("fixture $(entry.id) row normalization does not verify"))
    return true
end

function _assert_issue38_column_normalization(entry)
    normalizations = _require_field(entry, :normalizations)
    hasproperty(normalizations, :column) || throw(ArgumentError("fixture $(entry.id) missing column normalization metadata"))

    column = normalizations.column
    hasproperty(column, :core) || throw(ArgumentError("fixture $(entry.id) missing column normalization core"))
    hasproperty(column, :correction_factor) || throw(ArgumentError("fixture $(entry.id) missing column correction factor"))
    det(column.core) == one(base_ring(column.core)) || throw(ArgumentError("fixture $(entry.id) column normalization core must have determinant one"))
    entry.inputs.matrix * column.correction_factor == column.core || throw(ArgumentError("fixture $(entry.id) column normalization does not reconstruct the core"))
    return true
end

function _assert_issue38_expected_current_status(entry)
    status = _require_field(entry, :expected_current_status)

    if status isa AbstractString
        isempty(status) && throw(ArgumentError("fixture $(entry.id) expected current status must not be empty"))
        return true
    end

    if hasproperty(status, :status) || hasproperty(status, :state)
        current = hasproperty(status, :status) ? getproperty(status, :status) : getproperty(status, :state)
        current isa Symbol || throw(ArgumentError("fixture $(entry.id) expected current status must be symbolic"))
    end

    for field in (:message_substrings, :failure_substrings, :substrings, :messages)
        if hasproperty(status, field)
            substrings = getproperty(status, field)
            isempty(substrings) && throw(ArgumentError("fixture $(entry.id) expected current status must record failure substrings"))
            all(substring -> substring isa AbstractString, substrings) || throw(ArgumentError("fixture $(entry.id) expected current status substrings must be strings"))
            return true
        end
    end

    throw(ArgumentError("fixture $(entry.id) expected current status must include failure message substrings"))
end

function validate_toricbuilder_issue38_fixture(entry)
    _assert_issue38_metadata(entry)
    _assert_issue38_determinant_profile(entry)
    _assert_issue38_row_normalization(entry)
    _assert_issue38_column_normalization(entry)
    _assert_issue38_expected_current_status(entry)
    return true
end

function validate_toricbuilder_issue38_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("ToricBuilder Issue 38 fixture catalog missing cases"))
    isempty(catalog.cases) && throw(ArgumentError("ToricBuilder Issue 38 fixture catalog must not be empty"))
    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) || throw(ArgumentError("ToricBuilder Issue 38 fixture ids must be unique"))

    for entry in catalog.cases
        validate_toricbuilder_issue38_fixture(entry)
    end

    return true
end

@testset "ToricBuilder Issue 38 Q block fixture" begin
    @test isfile(TORICBUILDER_ISSUE38_FIXTURE_PATH)

    include(TORICBUILDER_ISSUE38_FIXTURE_PATH)
    catalog = ToricBuilderIssue38Cases.catalog()
    @test validate_toricbuilder_issue38_catalog(catalog)
end
