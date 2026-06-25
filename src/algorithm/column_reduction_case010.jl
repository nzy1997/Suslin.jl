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

function _laurent_unit_creation_candidate(column::Vector, R)
    _is_laurent_polynomial_ring(R) || return nothing
    length(column) == 5 || return nothing
    target_unit = one(R)

    for pivot_idx in eachindex(column), source_idx in eachindex(column)
        pivot_idx == source_idx && continue
        coeff = _try_laurent_divexact(target_unit - column[pivot_idx], column[source_idx])
        coeff === nothing && continue
        coeff == zero(R) && continue
        column[pivot_idx] + coeff * column[source_idx] == target_unit || continue
        return (;
            pivot_index = pivot_idx,
            source_index = source_idx,
            target_unit,
            creation_coefficient = coeff,
        )
    end

    return nothing
end

function _laurent_unit_creation_factors(n::Int, pivot_idx::Int, source_idx::Int, coeff, R)
    return [elementary_matrix(n, pivot_idx, source_idx, coeff, R)]
end

function _reduce_via_laurent_unit_creation_certificate(column::Vector, R)
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
    column::Vector,
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

function _reduce_laurent_unimodular_column_certificate(column::Vector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _unit_entry_reduction_certificate_stage(column, unit_idx, R)

    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return unit_creation

    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    is_unimodular_column(poly_column, P) || return nothing

    poly_result = _reduce_polynomial_unimodular_column_exact_certificate(poly_column, P)
    poly_result === nothing && return nothing

    polynomial_certificate = _ecp_certificate_from_stage(poly_column, P, poly_result.stage)
    lifted_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in polynomial_certificate.factors]
    shift = only(normalization.metadata.shift_monomials)
    inverse_shift = only(normalization.metadata.inverse_shift_monomials)
    normalization_factors = _unit_normalization_factors(length(column), inverse_shift, shift, R)
    factors = _checked_reduction_factors(
        vcat(normalization_factors, lifted_factors),
        column,
        R,
        "Laurent normalization reduction",
    )
    stage = (;
        kind = :laurent_normalization,
        input_column = _ecp_column_tuple(column),
        normalization,
        normalized_column = _ecp_column_tuple(poly_column),
        polynomial_certificate,
        lifted_factors,
        shift,
        inverse_shift,
        normalization_factors,
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

function _diagnose_laurent_unimodular_column_reduction(column::Vector, R, attempted::Vector{Symbol})
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return (; supported = true, stage = :unit_entry)

    push!(attempted, :laurent_unit_creation)
    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return (; supported = true, stage = :laurent_unit_creation)

    push!(attempted, :laurent_normalization)
    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    return _diagnose_polynomial_unimodular_column_reduction(poly_column, P, attempted)
end
