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

function _issue38_fixture_matrix_size(A)
    return (nrows(A), ncols(A))
end

function _require_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _assert_supported_column_peel_status(status, label::AbstractString)
    hasproperty(status, :kind) || throw(ArgumentError("fixture $(label) status missing kind"))
    status.kind == :supported_column_peel || throw(ArgumentError("fixture $(label) expected supported column peel status"))
    hasproperty(status, :factorization_path) || throw(ArgumentError("fixture $(label) status missing factorization path"))
    status.factorization_path == :laurent_column_peel || throw(ArgumentError("fixture $(label) expected Laurent column peel path"))
    return true
end

function _assert_supported_column_peel_core(matrix, status, label::AbstractString)
    _assert_supported_column_peel_status(status, label)
    factors = elementary_factorization(matrix)
    isempty(factors) && throw(ArgumentError("fixture $(label) expected nonempty factors"))
    verify_factorization(matrix, factors) || throw(ArgumentError("fixture $(label) factorization did not verify"))
    return true
end

function _assert_original_q_remains_unsupported(entry)
    err = try
        elementary_factorization(entry.inputs.matrix)
        nothing
    catch caught
        caught
    end

    err isa ArgumentError || throw(ArgumentError("fixture $(entry.id) original Q must remain unsupported"))
    message = sprint(showerror, err)
    occursin("elementary_factorization(A) is an elementary-only SL_n API", message) ||
        throw(ArgumentError("fixture $(entry.id) original Q determinant contract changed"))
    occursin("laurent_gl_factorization_certificate(A)", message) ||
        throw(ArgumentError("fixture $(entry.id) original Q certificate guidance changed"))
    return true
end

function _assert_issue38_metadata(entry)
    for field in REQUIRED_TORICBUILDER_ISSUE38_FIELDS
        _require_field(entry, field)
    end

    isempty(entry.consumer_test_ids) && throw(ArgumentError("fixture $(entry.id) must record at least one consumer test id"))
    entry.id == "toricbuilder-issue-38-q-block" || throw(ArgumentError("fixture $(entry.id) has unexpected id"))
    entry.kind == :toricbuilder_issue38_q_block || throw(ArgumentError("fixture $(entry.id) has unexpected kind $(entry.kind)"))
    hasproperty(entry.ring, :description) || throw(ArgumentError("fixture $(entry.id) missing ring description"))
    hasproperty(entry.ring, :object) || throw(ArgumentError("fixture $(entry.id) missing ring object"))
    hasproperty(entry.ring, :generators) || throw(ArgumentError("fixture $(entry.id) missing ring generators"))
    hasproperty(entry.ring, :variables) || throw(ArgumentError("fixture $(entry.id) missing ring variables"))
    hasproperty(entry.provenance, :issue) || throw(ArgumentError("fixture $(entry.id) missing provenance issue"))
    hasproperty(entry.provenance, :reported_main_commit) || throw(ArgumentError("fixture $(entry.id) missing provenance commit"))
    hasproperty(entry.inputs, :matrix) || throw(ArgumentError("fixture $(entry.id) missing matrix input"))
    hasproperty(entry.dimensions, :matrix) || throw(ArgumentError("fixture $(entry.id) missing matrix dimensions"))

    _issue38_fixture_matrix_size(entry.inputs.matrix) == (6, 6) || throw(ArgumentError("fixture $(entry.id) must be a 6x6 matrix"))
    entry.dimensions.matrix == (6, 6) || throw(ArgumentError("fixture $(entry.id) matrix dimensions must be 6x6"))

    u, v = entry.ring.generators
    profile = _require_field(entry, :determinant_profile)
    hasproperty(profile, :expected_class) || throw(ArgumentError("fixture $(entry.id) missing expected determinant class"))
    expected_determinant = _require_field(profile, :expected_determinant)
    expected_determinant == u * v || throw(ArgumentError("fixture $(entry.id) expected determinant must be u*v"))
    det(entry.inputs.matrix) == expected_determinant || throw(ArgumentError("fixture $(entry.id) determinant does not match metadata"))
    classify_laurent_determinant(entry.inputs.matrix).classification == profile.expected_class || throw(ArgumentError("fixture $(entry.id) determinant class does not match metadata"))
    return true
end

function _assert_issue38_determinant_profile(entry)
    profile = _require_field(entry, :determinant_profile)
    hasproperty(profile, :expected_determinant) || throw(ArgumentError("fixture $(entry.id) missing expected determinant"))
    hasproperty(profile, :expected_class) || throw(ArgumentError("fixture $(entry.id) missing expected determinant class"))
    return true
end

function _assert_issue38_row_normalization(entry)
    normalizations = _require_field(entry, :normalizations)
    hasproperty(normalizations, :row) || throw(ArgumentError("fixture $(entry.id) missing row normalization metadata"))

    row = normalizations.row
    hasproperty(row, :core) || throw(ArgumentError("fixture $(entry.id) missing row normalization core"))
    hasproperty(row, :normalization) || throw(ArgumentError("fixture $(entry.id) missing row normalization object"))
    hasproperty(row, :expected_current_status) || throw(ArgumentError("fixture $(entry.id) missing row expected current status"))
    row.core == row.normalization.normalized_matrix || throw(ArgumentError("fixture $(entry.id) row core must equal normalized matrix"))
    det(row.core) == one(base_ring(row.core)) || throw(ArgumentError("fixture $(entry.id) row normalization core must have determinant one"))
    verify_laurent_gl_normalization(entry.inputs.matrix, row.normalization) || throw(ArgumentError("fixture $(entry.id) row normalization does not verify"))
    _assert_supported_column_peel_core(row.normalization.normalized_matrix, row.expected_current_status, "$(entry.id) row")
    return true
end

function _assert_issue38_column_normalization(entry)
    normalizations = _require_field(entry, :normalizations)
    hasproperty(normalizations, :column) || throw(ArgumentError("fixture $(entry.id) missing column normalization metadata"))

    column = normalizations.column
    hasproperty(column, :core) || throw(ArgumentError("fixture $(entry.id) missing column normalization core"))
    hasproperty(column, :correction_factor) || throw(ArgumentError("fixture $(entry.id) missing column correction factor"))
    hasproperty(column, :expected_current_status) || throw(ArgumentError("fixture $(entry.id) missing column expected current status"))
    det(column.core) == one(base_ring(column.core)) || throw(ArgumentError("fixture $(entry.id) column normalization core must have determinant one"))
    entry.inputs.matrix * column.correction_factor == column.core || throw(ArgumentError("fixture $(entry.id) column normalization does not reconstruct the core"))
    _assert_supported_column_peel_core(column.core, column.expected_current_status, "$(entry.id) column")
    return true
end

function _assert_issue38_expected_current_status(entry)
    status = _require_field(entry, :expected_current_status)
    hasproperty(status, :row) || throw(ArgumentError("fixture $(entry.id) expected current status missing row metadata"))
    hasproperty(status, :column) || throw(ArgumentError("fixture $(entry.id) expected current status missing column metadata"))
    _assert_supported_column_peel_status(status.row, "$(entry.id) expected_current_status.row")
    _assert_supported_column_peel_status(status.column, "$(entry.id) expected_current_status.column")
    return true
end

function validate_toricbuilder_issue38_fixture(entry)
    _assert_issue38_metadata(entry)
    _assert_issue38_determinant_profile(entry)
    _assert_issue38_row_normalization(entry)
    _assert_issue38_column_normalization(entry)
    _assert_issue38_expected_current_status(entry)
    _assert_original_q_remains_unsupported(entry)
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

    entry = only(catalog.cases)
    corrupted_matrix = copy(entry.inputs.matrix)
    corrupted_matrix[1, 1] += one(base_ring(corrupted_matrix))
    bad_entry = merge(entry, (; inputs = merge(entry.inputs, (; matrix = corrupted_matrix))))
    @test_throws ArgumentError validate_toricbuilder_issue38_fixture(bad_entry)
end
