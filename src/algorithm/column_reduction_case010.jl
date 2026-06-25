function _try_laurent_divexact(numerator, denominator)
    iszero(denominator) && return nothing
    try
        return divexact(numerator, denominator)
    catch err
        err isa InterruptException && rethrow()
        (err isa ErrorException || err isa ArgumentError || err isa MethodError) && return nothing
        rethrow()
    end
end

function _preferred_laurent_unit_creation_candidate(left, right)
    left === nothing && return right
    # Prefer lower rows so the recursive column-peel keeps the exposed case_010
    # trailing blocks in the existing supported Laurent families.
    right.pivot_index > left.pivot_index && return right
    right.pivot_index == left.pivot_index && right.source_index < left.source_index && return right
    return left
end

function _laurent_unit_creation_candidate(column::AbstractVector, R)
    _is_laurent_polynomial_ring(R) || return nothing
    length(column) in (4, 5) || return nothing
    target_unit = one(R)
    candidate = nothing

    for pivot_idx in eachindex(column), source_idx in eachindex(column)
        pivot_idx == source_idx && continue
        coeff = _try_laurent_divexact(target_unit - column[pivot_idx], column[source_idx])
        coeff === nothing && continue
        coeff == zero(R) && continue
        column[pivot_idx] + coeff * column[source_idx] == target_unit || continue
        candidate = _preferred_laurent_unit_creation_candidate(candidate, (;
            pivot_index = pivot_idx,
            source_index = source_idx,
            target_unit,
            creation_coefficient = coeff,
        ))
    end

    return candidate
end

function _laurent_unit_creation_factors(n::Int, pivot_idx::Int, source_idx::Int, coeff, R)
    return [elementary_matrix(n, pivot_idx, source_idx, coeff, R)]
end

function _reduce_via_laurent_unit_creation_certificate(column::AbstractVector, R)
    candidate = _laurent_unit_creation_candidate(column, R)
    candidate === nothing && return nothing
    return _laurent_unit_creation_certificate_stage(
        column,
        R,
        candidate.pivot_index,
        candidate.source_index,
        candidate.target_unit,
        candidate.creation_coefficient,
    )
end

function _laurent_unit_creation_certificate_stage(
    column::AbstractVector,
    R,
    pivot_idx::Int,
    source_idx::Int,
    target_unit,
    creation_coefficient,
)
    n = length(column)
    creation_factors = _laurent_unit_creation_factors(n, pivot_idx, source_idx, creation_coefficient, R)
    created_column_matrix = _apply_reduction_factors(creation_factors, column, R)
    created_column = collect(_ecp_matrix_column_to_tuple(created_column_matrix))
    created_column[pivot_idx] == target_unit ||
        error("internal Laurent unit-creation stage did not create the target unit")

    unit_stage = _unit_entry_reduction_certificate_stage(created_column, pivot_idx, R).stage
    factors = _checked_reduction_factors(
        vcat(unit_stage.factors, creation_factors),
        column,
        R,
        "Laurent unit-creation reduction",
    )
    stage = (;
        kind = :laurent_unit_creation,
        input_column = _ecp_column_tuple(column),
        pivot_index = pivot_idx,
        source_index = source_idx,
        target_unit,
        creation_coefficient,
        creation_factors,
        created_column = tuple(created_column...),
        unit_stage,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end

function _ecp_replay_stage(
    stage::NamedTuple{
        (
            :kind,
            :input_column,
            :pivot_index,
            :source_index,
            :target_unit,
            :creation_coefficient,
            :creation_factors,
            :created_column,
            :unit_stage,
            :factors,
            :output_column,
        ),
        T,
    },
    input_column,
    R,
) where {T}
    invalid_replay = (; ok = false, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
    stage.kind == :laurent_unit_creation || return invalid_replay
    _is_laurent_polynomial_ring(R) || return invalid_replay

    pivot_idx = stage.pivot_index
    source_idx = stage.source_index
    pivot_idx isa Integer && source_idx isa Integer || return invalid_replay
    1 <= pivot_idx <= length(input_column) || return invalid_replay
    1 <= source_idx <= length(input_column) || return invalid_replay
    pivot_idx == source_idx && return invalid_replay
    stage.target_unit == one(R) || return invalid_replay

    coeff = _try_laurent_divexact(stage.target_unit - input_column[pivot_idx], input_column[source_idx])
    coeff === nothing && return invalid_replay
    creation_factors = _laurent_unit_creation_factors(length(input_column), pivot_idx, source_idx, coeff, R)
    created_column = collect(_ecp_matrix_column_to_tuple(_apply_reduction_factors(creation_factors, input_column, R)))
    unit_replay = _ecp_replay_stage(stage.unit_stage, created_column, R)
    expected_factors = vcat(unit_replay.factors, creation_factors)
    expected_output = _apply_reduction_factors(expected_factors, input_column, R)
    ok = stage.input_column == _ecp_column_tuple(input_column) &&
        stage.creation_coefficient == coeff &&
        created_column[pivot_idx] == stage.target_unit &&
        stage.created_column == tuple(created_column...) &&
        unit_replay.ok &&
        _ecp_factor_sequences_equal(stage.creation_factors, creation_factors) &&
        _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
        stage.output_column == expected_output
    return (; ok, factors = expected_factors, output_column = expected_output)
end
