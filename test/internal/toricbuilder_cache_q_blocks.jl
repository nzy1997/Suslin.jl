using Test
using Suslin
using Oscar

const TORICBUILDER_CACHE_Q_BLOCKS_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_cache_q_blocks.jl")
const TORICBUILDER_ISSUE38_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_issue38_cases.jl")

const REQUIRED_TORICBUILDER_CACHE_Q_BLOCK_FIELDS = (
    :id,
    :kind,
    :source_cache_file,
    :source_block,
    :ring,
    :dimensions,
    :storage_format,
    :value_format,
    :sparse_entry_count,
    :sparse_entries,
    :matrix,
    :expected_test_level,
    :expected_current_status,
    :provenance,
)

function _cache_qblock_field(entry, field::Symbol)
    hasproperty(entry, field) ||
        throw(ArgumentError("ToricBuilder cache Q-block entry $(entry.id) missing field $(field)"))
    return getproperty(entry, field)
end

function _cache_qblock_matrix_size(A)
    return (nrows(A), ncols(A))
end

function validate_toricbuilder_cache_qblock_entry(entry)
    for field in REQUIRED_TORICBUILDER_CACHE_Q_BLOCK_FIELDS
        _cache_qblock_field(entry, field)
    end

    startswith(entry.id, "case_") ||
        throw(ArgumentError("ToricBuilder cache Q-block id must start with case_"))
    entry.kind == :toricbuilder_cache_column_q_block ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) has unexpected kind $(entry.kind)"))
    entry.source_block == :column_transformation_upper_left_q_block ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) has unexpected source block"))
    entry.storage_format == :sparse_coordinate_entries ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) must use sparse coordinate storage"))
    entry.value_format == :julia_laurent_expression_string ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) has unexpected value format"))
    entry.sparse_entry_count == length(entry.sparse_entries) ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) sparse entry count mismatch"))
    entry.sparse_entry_count > 0 ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) must not be empty"))
    hasproperty(entry.dimensions, :matrix) ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) missing matrix dimensions"))
    hasproperty(entry.dimensions, :column_transformation) ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) missing source dimensions"))
    entry.dimensions.column_transformation[1] == 2 * entry.dimensions.matrix[1] ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) source row dimension mismatch"))
    entry.dimensions.column_transformation[2] == 2 * entry.dimensions.matrix[2] ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) source column dimension mismatch"))
    hasproperty(entry.ring, :variables) ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) missing ring variables"))
    entry.ring.variables == ("u", "v") ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) expected u/v variables"))
    hasproperty(entry.matrix, :storage) ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) missing matrix storage descriptor"))
    entry.matrix.storage == :lazy_sparse_coordinate_entries ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) must lazily materialize the matrix"))
    entry.matrix.materializer == :materialize_matrix ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) has unexpected matrix materializer"))

    seen = Set{Tuple{Int, Int}}()
    for (row, col, value) in entry.sparse_entries
        1 <= row <= entry.dimensions.matrix[1] ||
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) sparse row out of range"))
        1 <= col <= entry.dimensions.matrix[2] ||
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) sparse column out of range"))
        value isa AbstractString ||
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) sparse value must be a string"))
        !isempty(value) ||
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) sparse value must not be empty"))
        !occursin("AbstractAlgebra.Generic", value) ||
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) sparse value must not use internal constructors"))
        (row, col) in seen &&
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) duplicate sparse coordinate"))
        push!(seen, (row, col))
    end

    if entry.dimensions.matrix[1] >= 100
        entry.expected_test_level == :optional_slow ||
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) large matrix must be optional slow"))
    end

    return true
end

function validate_toricbuilder_cache_qblock_materialization(fixture_module, entry)
    matrix = fixture_module.materialize_matrix(entry)
    _cache_qblock_matrix_size(matrix) == entry.dimensions.matrix ||
        throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) materialized matrix dimensions mismatch"))

    R = base_ring(matrix)
    for (row, col, _) in entry.sparse_entries
        matrix[row, col] != zero(R) ||
            throw(ArgumentError("ToricBuilder cache Q-block $(entry.id) sparse entry materialized to zero"))
    end

    return matrix
end

function validate_toricbuilder_cache_qblock_catalog(catalog)
    hasproperty(catalog, :cases) ||
        throw(ArgumentError("ToricBuilder cache Q-block catalog missing cases"))
    length(catalog.cases) >= 12 ||
        throw(ArgumentError("ToricBuilder cache Q-block catalog must include at least case_001 through case_012"))

    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) ||
        throw(ArgumentError("ToricBuilder cache Q-block ids must be unique"))
    required_ids = Set(["case_$(lpad(string(idx), 3, "0"))" for idx in 1:12])
    required_ids ⊆ Set(ids) ||
        throw(ArgumentError("ToricBuilder cache Q-block catalog missing required case IDs"))

    for entry in catalog.cases
        validate_toricbuilder_cache_qblock_entry(entry)
    end

    return true
end

@testset "ToricBuilder cache Q-block fixtures" begin
    @test isfile(TORICBUILDER_CACHE_Q_BLOCKS_PATH)

    include(TORICBUILDER_CACHE_Q_BLOCKS_PATH)
    catalog = ToricBuilderCacheQBlocks.catalog()
    @test validate_toricbuilder_cache_qblock_catalog(catalog)

    by_id = Dict(entry.id => entry for entry in catalog.cases)
    @test by_id["case_001"].dimensions.matrix == (6, 6)
    @test by_id["case_011"].dimensions.matrix == (288, 288)
    @test by_id["case_011"].expected_test_level == :optional_slow

    case001_matrix =
        validate_toricbuilder_cache_qblock_materialization(ToricBuilderCacheQBlocks, by_id["case_001"])

    include(TORICBUILDER_ISSUE38_FIXTURE_PATH)
    issue38 = only(ToricBuilderIssue38Cases.catalog().cases)
    @test case001_matrix == issue38.inputs.matrix

    bad_dimensions = merge(
        by_id["case_001"],
        (; dimensions = merge(by_id["case_001"].dimensions, (; matrix = (5, 5),))),
    )
    @test_throws ArgumentError validate_toricbuilder_cache_qblock_entry(bad_dimensions)
end
