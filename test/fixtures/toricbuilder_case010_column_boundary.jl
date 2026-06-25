module ToricBuilderCase010ColumnBoundary

using Oscar
using Suslin

include(joinpath(@__DIR__, "toricbuilder_cache_q_blocks.jl"))

const CASE_ID = "case_010"
const EXPECTED_RING_DESCRIPTION = "GF(2)[u^+/-1, v^+/-1]"
const EXPECTED_DIAGNOSTIC = "unsupported exact unimodular column reduction"
const REQUIRED_BOUNDARY_FIELDS = (
    :case_id,
    :original_matrix,
    :normalization,
    :normalized_matrix,
    :first_failing_peel_dimension,
    :failing_input_matrix,
    :failing_column,
    :ring,
    :ring_description,
    :expected_diagnostic,
)

function _case010_entry()
    matches = filter(entry -> entry.id == CASE_ID, ToricBuilderCacheQBlocks.catalog().cases)
    length(matches) == 1 ||
        throw(ArgumentError("expected exactly one ToricBuilder cache Q-block entry for $(CASE_ID)"))
    return only(matches)
end

function _ring_generator_names(R)
    return Tuple(string.(gens(R)))
end

function _is_expected_uv_laurent_ring(R)::Bool
    try
        return Suslin._is_laurent_polynomial_ring(R) &&
            base_ring(R) == GF(2) &&
            _ring_generator_names(R) == ("u", "v")
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _column_entries_are_over_ring(column, R)::Bool
    try
        return all(entry -> R(entry) == entry, column)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _first_failing_boundary(normalized_matrix)
    nrows(normalized_matrix) == 6 ||
        throw(ArgumentError("$(CASE_ID) expected a normalized 6x6 matrix before extracting the preserved 5x5 boundary"))
    step = Suslin._laurent_column_peel_step(normalized_matrix)
    current = step.next_block
    d = nrows(current)
    d == 5 || throw(ArgumentError("$(CASE_ID) expected the preserved boundary block to have dimension 5"))
    return (;
        first_failing_peel_dimension = d,
        failing_input_matrix = current,
        failing_column = [current[row, d] for row in 1:d],
        passed_peel_dimensions = (6,),
        observed_diagnostic = nothing,
    )
end

function boundary_fixture()
    entry = _case010_entry()
    original_matrix = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    normalization = Suslin.normalize_laurent_gl_matrix(original_matrix)
    boundary = _first_failing_boundary(normalization.normalized_matrix)

    return (;
        case_id = entry.id,
        source_entry = entry,
        original_matrix,
        normalization,
        normalized_matrix = normalization.normalized_matrix,
        ring = base_ring(original_matrix),
        ring_description = entry.ring.description,
        expected_diagnostic = EXPECTED_DIAGNOSTIC,
        boundary...,
    )
end

function _has_required_boundary_fields(fixture)::Bool
    return all(field -> hasproperty(fixture, field), REQUIRED_BOUNDARY_FIELDS)
end

function validate_boundary_fixture(fixture)::Symbol
    _has_required_boundary_fields(fixture) || return :missing_metadata
    fixture.case_id == CASE_ID || return :wrong_case
    fixture.first_failing_peel_dimension == 5 || return :wrong_peel_dimension
    length(fixture.failing_column) == 5 || return :wrong_column_length
    fixture.ring_description == EXPECTED_RING_DESCRIPTION || return :wrong_ring

    R = fixture.ring
    _is_expected_uv_laurent_ring(R) || return :wrong_ring
    _column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
    Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular

    return :ok
end

function non_unimodular_negative_control(fixture = boundary_fixture())
    _, v = gens(fixture.ring)
    nonunit = v + one(fixture.ring)
    return merge(fixture, (; failing_column = [nonunit * entry for entry in fixture.failing_column]))
end

function single_entry_zero_perturbations(fixture = boundary_fixture())
    return [
        begin
            column = copy(fixture.failing_column)
            column[idx] = zero(fixture.ring)
            column
        end
        for idx in eachindex(fixture.failing_column)
    ]
end

end
