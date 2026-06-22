struct ECPColumnReductionCertificate
    original_column
    ring
    stages
    factors::Vector
    final_column
    verification
end

struct ECPMonicitySearchResult
    original_column
    ring
    variable_order
    max_shift_power::Int
    shift_signs
    source_variable_index::Int
    source_variable
    target_variable_index::Int
    target_variable
    shift_power::Int
    shift_sign
    shift_polynomial
    selected_monic_index::Int
    selected_monic_entry
    stage
    factors::Vector
end

struct ECPMonicitySearchFailure
    kind::Symbol
    original_column
    ring
    variable_order
    max_shift_power::Int
    shift_signs
    source_variables
    target_variable
    shift_powers
    shift_polynomials
    attempted_candidates::Int
    message::String
end

function reduce_unimodular_column(v::AbstractVector, R)
    return ecp_column_reduction_certificate(v, R).factors
end

function ecp_column_reduction_certificate(v::AbstractVector, R)
    column = _validated_unimodular_column(v, R)
    result = _is_laurent_polynomial_ring(R) ?
        _reduce_laurent_unimodular_column_certificate(column, R) :
        _reduce_polynomial_unimodular_column_exact_certificate(column, R)
    result !== nothing || _throw_unsupported_unimodular_column_reduction(column, R)

    factors = _checked_reduction_factors(result.factors, column, R, "certificate reducer")
    stages = ((; kind = :validation, input_length = length(column), is_unimodular = true), result.stage)
    final_column = _apply_reduction_factors(factors, column, R)
    provisional = ECPColumnReductionCertificate(column, R, stages, factors, final_column, nothing)
    verification = _ecp_column_reduction_replay_summary(provisional)
    verification.overall_ok || error("internal ECP column reduction certificate verification failed")
    certificate = ECPColumnReductionCertificate(column, R, stages, factors, final_column, verification)
    verify_ecp_column_reduction(certificate) || error("internal ECP column reduction certificate storage verification failed")
    return certificate
end

function _validated_unimodular_column(v::AbstractVector, R)
    Base.require_one_based_indexing(v)
    n = length(v)
    n >= 3 || throw(ArgumentError("v must have length at least 3"))

    column = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    is_unimodular_column(column, R) || throw(ArgumentError("v must be a unimodular column"))
    return column
end

function _factor_sequence_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        nrows(factor) == n || throw(ArgumentError("reduction factor has wrong row count"))
        ncols(factor) == n || throw(ArgumentError("reduction factor has wrong column count"))
        _same_base_ring(base_ring(factor), R) || throw(ArgumentError("reduction factor has wrong base ring"))
        product *= factor
    end
    return product
end

function _apply_reduction_factors(factors, column::AbstractVector, R)
    n = length(column)
    return _factor_sequence_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _target_reduced_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _checked_reduction_factors(factors, column::AbstractVector, R, stage::AbstractString)
    n = length(column)
    _apply_reduction_factors(factors, column, R) == _target_reduced_column(R, n) ||
        throw(ErrorException("internal error: $(stage) produced factors that do not reduce the column to e_n"))
    return factors
end

_reduce_unit_witness_column(column::AbstractVector, pivot_idx::Int, R) = _unit_entry_reduction_certificate_stage(column, pivot_idx, R).factors

function _unit_entry_reduction_certificate_stage(column::AbstractVector, pivot_idx::Int, R)
    n = length(column)
    u = column[pivot_idx]
    uinv = inv(u)
    elimination_factors = typeof(identity_matrix(R, n))[]

    for k in 1:n
        k == pivot_idx && continue
        coeff = -column[k] * uinv
        coeff == zero(R) && continue
        push!(elimination_factors, elementary_matrix(n, k, pivot_idx, coeff, R))
    end

    move_factors = typeof(identity_matrix(R, n))[]
    if pivot_idx != n
        push!(move_factors, elementary_matrix(n, pivot_idx, n, -one(R), R))
        push!(move_factors, elementary_matrix(n, n, pivot_idx, one(R), R))
    end

    normalization_factors = _unit_normalization_factors(n, u, uinv, R)
    factors = _checked_reduction_factors(
        vcat(normalization_factors, move_factors, elimination_factors),
        column,
        R,
        "unit entry reduction",
    )
    stage = (;
        kind = :unit_entry,
        input_column = _ecp_column_tuple(column),
        pivot_index = pivot_idx,
        pivot_value = column[pivot_idx],
        pivot_inverse = uinv,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end

function _unimodular_witness(column::AbstractVector, R)
    witness_matrix = Oscar.coordinates(one(R), ideal(R, column))
    return [witness_matrix[1, idx] for idx in 1:ncols(witness_matrix)]
end

function _witness_unit_creation_factors(column::AbstractVector, witness::AbstractVector, pivot_idx::Int, R)
    n = length(column)
    u = witness[pivot_idx]
    uinv = inv(u)
    unit_creation_factors = typeof(identity_matrix(R, n))[]

    for k in 1:n
        k == pivot_idx && continue
        coeff = uinv * witness[k]
        coeff == zero(R) && continue
        push!(unit_creation_factors, elementary_matrix(n, pivot_idx, k, coeff, R))
    end

    return unit_creation_factors
end

function _reduce_via_witness_unit(column::AbstractVector, witness::AbstractVector, pivot_idx::Int, R)
    return _witness_unit_reduction_certificate_stage(column, witness, pivot_idx, R).factors
end

function _witness_unit_reduction_certificate_stage(column::AbstractVector, witness::AbstractVector, pivot_idx::Int, R)
    unit_creation_factors = _witness_unit_creation_factors(column, witness, pivot_idx, R)
    created_column_matrix = _apply_reduction_factors(unit_creation_factors, column, R)
    created_column = collect(_ecp_matrix_column_to_tuple(created_column_matrix))
    unit_stage = _unit_entry_reduction_certificate_stage(created_column, pivot_idx, R).stage
    factors = _checked_reduction_factors(
        vcat(unit_stage.factors, unit_creation_factors),
        column,
        R,
        "witness unit reduction",
    )
    stage = (;
        kind = :witness_unit,
        input_column = _ecp_column_tuple(column),
        witness = _ecp_column_tuple(witness),
        pivot_index = pivot_idx,
        witness_unit = witness[pivot_idx],
        witness_unit_inverse = inv(witness[pivot_idx]),
        unit_creation_factors,
        created_column = tuple(created_column...),
        unit_stage,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end

function _reduce_supported_unimodular_column(column::AbstractVector, R)
    result = _reduce_supported_unimodular_column_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_supported_unimodular_column_certificate(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _unit_entry_reduction_certificate_stage(column, unit_idx, R)

    witness = _unimodular_witness(column, R)
    witness_unit_idx = findfirst(is_unit, witness)
    witness_unit_idx !== nothing && return _witness_unit_reduction_certificate_stage(column, witness, witness_unit_idx, R)

    return nothing
end

function _reduce_polynomial_unimodular_column_exact(column::AbstractVector, R)
    result = _reduce_polynomial_unimodular_column_exact_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_polynomial_unimodular_column_exact_certificate(column::AbstractVector, R)
    factors = _reduce_exact_small_column_certificate(column, R)
    factors !== nothing && return factors

    if length(column) > 3
        block_factors = _reduce_via_supported_three_block_certificate(column, R)
        block_factors !== nothing && return block_factors
    end

    return nothing
end

function _has_at_least_two_generators(R)::Bool
    try
        return ngens(R) >= 2
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        return false
    end
end

function _reduce_exact_small_column(column::AbstractVector, R)
    result = _reduce_exact_small_column_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_exact_small_column_certificate(column::AbstractVector, R)
    factors = _reduce_supported_unimodular_column_certificate(column, R)
    factors !== nothing && return factors

    _has_at_least_two_generators(R) || return nothing
    normalization_factors = _reduce_after_monicity_normalization_certificate(column, R)
    normalization_factors !== nothing && return normalization_factors

    return nothing
end

function _reduce_via_supported_three_block(column::AbstractVector, R)
    result = _reduce_via_supported_three_block_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_via_supported_three_block_certificate(column::AbstractVector, R)
    n = length(column)
    n > 3 || return nothing

    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        indices = (i, j, k)
        subcolumn = [column[idx] for idx in indices]
        is_unimodular_column(subcolumn, R) || continue

        subresult = _reduce_exact_small_column_certificate(subcolumn, R)
        subresult === nothing && continue

        return _embedded_three_block_reduction_certificate_stage(column, R, indices, subresult)
    end

    return nothing
end

function _embedded_three_block_reduction(column::AbstractVector, R, indices, subfactors)
    n = length(column)
    pivot_idx = indices[end]
    embedded_factors = [block_embedding(factor, n, indices) for factor in subfactors]
    after_block = _apply_reduction_factors(embedded_factors, column, R)

    elimination_factors = typeof(identity_matrix(R, n))[]
    for row in 1:n
        row == pivot_idx && continue
        coeff = -after_block[row, 1]
        coeff == zero(R) && continue
        push!(elimination_factors, elementary_matrix(n, row, pivot_idx, coeff, R))
    end

    move_factors = typeof(identity_matrix(R, n))[]
    if pivot_idx != n
        push!(move_factors, elementary_matrix(n, pivot_idx, n, -one(R), R))
        push!(move_factors, elementary_matrix(n, n, pivot_idx, one(R), R))
    end

    return _checked_reduction_factors(
        vcat(move_factors, elimination_factors, embedded_factors),
        column,
        R,
        "embedded 3-entry reduction",
    )
end

function _embedded_three_block_reduction_certificate_stage(column::AbstractVector, R, indices, subresult)
    n = length(column)
    subcolumn = [column[idx] for idx in indices]
    subcertificate = _ecp_certificate_from_stage(subcolumn, R, subresult.stage)
    pivot_idx = indices[end]
    embedded_factors = [block_embedding(factor, n, indices) for factor in subcertificate.factors]
    after_block = _apply_reduction_factors(embedded_factors, column, R)

    elimination_factors = typeof(identity_matrix(R, n))[]
    for row in 1:n
        row == pivot_idx && continue
        coeff = -after_block[row, 1]
        coeff == zero(R) && continue
        push!(elimination_factors, elementary_matrix(n, row, pivot_idx, coeff, R))
    end

    move_factors = typeof(identity_matrix(R, n))[]
    if pivot_idx != n
        push!(move_factors, elementary_matrix(n, pivot_idx, n, -one(R), R))
        push!(move_factors, elementary_matrix(n, n, pivot_idx, one(R), R))
    end

    factors = _checked_reduction_factors(
        vcat(move_factors, elimination_factors, embedded_factors),
        column,
        R,
        "embedded 3-entry reduction",
    )
    stage = (;
        kind = :embedded_three_block,
        input_column = _ecp_column_tuple(column),
        indices = tuple(indices...),
        subcolumn = _ecp_column_tuple(subcolumn),
        subcertificate,
        embedded_factors,
        post_block_column = _ecp_matrix_column_to_tuple(after_block),
        elimination_factors,
        move_factors,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end

function _reduce_after_monicity_normalization(column::AbstractVector, R)
    result = _reduce_after_monicity_normalization_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_after_monicity_normalization_certificate(column::AbstractVector, R)
    result = _deterministic_ecp_monicity_search(column, R)
    return result isa ECPMonicitySearchResult ? (; factors = result.factors, stage = result.stage) : nothing
end

function _deterministic_ecp_monicity_search(
    column::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    max_shift_power::Integer = 3,
    shift_signs = (one(R), -one(R)),
)
    normalized_order = _ecp_normalize_variable_order(R, variable_order)
    max_shift_power < 0 && throw(ArgumentError("max_shift_power must be nonnegative"))

    signs = tuple((_coerce_into_ring(R, sign, "shift sign") for sign in shift_signs)...)

    isempty(normalized_order) && return _ecp_monicity_search_failure(column, R, normalized_order, max_shift_power, signs, 0)
    length(normalized_order) < 2 && return _ecp_monicity_search_failure(column, R, normalized_order, max_shift_power, signs, 0)

    target_variable = normalized_order[end]
    target_variable_index = _ecp_generator_index(R, target_variable)
    attempted = 0
    for source_variable in normalized_order[1:(end - 1)]
        source_variable_index = _ecp_generator_index(R, source_variable)
        for shift_power in 1:max_shift_power, shift_sign in signs
            attempted += 1
            candidate = _ecp_monicity_candidate_stage(
                column,
                R,
                tuple(normalized_order...),
                source_variable_index,
                target_variable_index,
                shift_power,
                shift_sign,
            )
            candidate === nothing && continue
            stage = candidate.stage
            return ECPMonicitySearchResult(
                tuple(column...),
                R,
                tuple(normalized_order...),
                Int(max_shift_power),
                signs,
                stage.source_variable_index,
                stage.source_variable,
                stage.target_variable_index,
                stage.target_variable,
                stage.shift_power,
                stage.shift_sign,
                stage.shift_polynomial,
                stage.selected_monic_index,
                stage.selected_monic_entry,
                stage,
                candidate.factors,
            )
        end
    end

    return _ecp_monicity_search_failure(column, R, normalized_order, max_shift_power, signs, attempted)
end

function _ecp_normalize_variable_order(R, variable_order)
    ring_gens = collect(gens(R))
    normalized = Any[]
    for variable in variable_order
        match_idx = _ecp_variable_order_match_index(ring_gens, variable)
        match_idx === nothing && throw(ArgumentError("variable_order must contain only generators of R"))
        generator = ring_gens[match_idx]
        any(existing -> existing == generator, normalized) &&
            throw(ArgumentError("variable_order must not contain duplicate generators"))
        push!(normalized, generator)
    end
    return tuple(normalized...)
end

function _ecp_variable_order_match_index(ring_gens::AbstractVector, variable)
    match_idx = findfirst(generator -> generator == variable, ring_gens)
    match_idx !== nothing && return match_idx
    variable isa Symbol || return nothing
    return findfirst(generator -> Symbol(string(generator)) == variable, ring_gens)
end

function _ecp_generator_index(R, variable)
    generator_idx = findfirst(generator -> generator == variable, gens(R))
    generator_idx === nothing && throw(ArgumentError("variable order must contain only generators of R"))
    return generator_idx
end

function _ecp_monicity_candidate_stage(
    column::AbstractVector,
    R,
    variable_order,
    source_variable_index::Int,
    target_variable_index::Int,
    shift_power::Int,
    shift_sign,
)
    ring_gens = collect(gens(R))
    source_variable = ring_gens[source_variable_index]
    target_variable = ring_gens[target_variable_index]

    forward_values = copy(ring_gens)
    forward_values[source_variable_index] = source_variable + shift_sign * target_variable^shift_power
    transformed = [_coerce_into_ring(R, evaluate(entry, forward_values), "substituted column entry") for entry in column]
    selected_monic_index = _ecp_first_monic_entry_index(transformed, R, target_variable_index)
    selected_monic_index === nothing && return nothing

    transformed_result = _reduce_supported_unimodular_column_certificate(transformed, R)
    transformed_result === nothing && return nothing

    inverse_values = copy(ring_gens)
    inverse_values[source_variable_index] = source_variable - shift_sign * target_variable^shift_power
    inverse_substituted_factors = [
        _substitute_matrix_entries(factor, inverse_values, R)
        for factor in transformed_result.factors
    ]
    factors = _checked_reduction_factors(
        inverse_substituted_factors,
        column,
        R,
        "monicity normalization reduction",
    )
    target = _target_reduced_column(R, length(column))
    transformed_output = _apply_reduction_factors(transformed_result.factors, transformed, R)
    original_output = _apply_reduction_factors(factors, column, R)
    forward_substitution = _ecp_substitution_map_tuple(ring_gens, forward_values)
    inverse_substitution = _ecp_substitution_map_tuple(ring_gens, inverse_values)
    selected_monic_entry = transformed[selected_monic_index]
    first_coordinate_move_factors = typeof(identity_matrix(R, length(column)))[]
    first_coordinate_column = tuple(transformed...)
    variable_change_verification = (;
        selected_monic_ok = _is_monic_in_variable(selected_monic_entry, R, target_variable_index),
        transformed_reduction_ok = transformed_output == target,
        original_reduction_ok = original_output == target,
    )
    stage = (;
        kind = :monicity_normalization,
        input_column = _ecp_column_tuple(column),
        variable_order = tuple(variable_order...),
        variable_index = source_variable_index,
        source_variable_index,
        source_variable,
        last_variable_index = target_variable_index,
        target_variable_index,
        target_variable,
        shift_power,
        shift_sign,
        shift_polynomial = shift_sign * target_variable^shift_power,
        forward_values = tuple(forward_values...),
        inverse_values = tuple(inverse_values...),
        forward_substitution,
        inverse_substitution,
        transformed_column = tuple(transformed...),
        selected_monic_index,
        selected_monic_entry,
        first_coordinate_strategy = _ecp_first_coordinate_strategy(selected_monic_index),
        first_coordinate_move_factors,
        first_coordinate_column,
        transformed_stage = transformed_result.stage,
        transformed_factors = transformed_result.factors,
        inverse_substituted_factors,
        factors,
        output_column = original_output,
        variable_change_verification,
    )
    return (; factors, stage)
end

function _ecp_monicity_search_failure(column::AbstractVector, R, variable_order, max_shift_power::Integer, shift_signs, attempted_candidates::Int)
    source_variables = length(variable_order) < 2 ? () : tuple(variable_order[1:(end - 1)]...)
    target_variable = isempty(variable_order) ? nothing : variable_order[end]
    shift_powers = tuple((1:Int(max_shift_power))...)
    shift_polynomials = target_variable === nothing ? () : tuple((
        sign * target_variable^power
        for power in shift_powers
        for sign in shift_signs
    )...)
    message = "exhausted deterministic ECP monicity search for variable_order=$(tuple(variable_order...)), max_shift_power=$(Int(max_shift_power)), attempted_candidates=$(attempted_candidates)"
    return ECPMonicitySearchFailure(
        :monicity_search_exhausted,
        tuple(column...),
        R,
        tuple(variable_order...),
        Int(max_shift_power),
        tuple(shift_signs...),
        source_variables,
        target_variable,
        shift_powers,
        shift_polynomials,
        attempted_candidates,
        message,
    )
end

function _reduce_laurent_unimodular_column(column::AbstractVector, R)
    result = _reduce_laurent_unimodular_column_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_laurent_unimodular_column_certificate(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _unit_entry_reduction_certificate_stage(column, unit_idx, R)

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

function _ecp_certificate_from_stage(column::AbstractVector, R, stage)
    factors = _checked_reduction_factors(stage.factors, column, R, "certificate stage")
    stages = ((; kind = :validation, input_length = length(column), is_unimodular = true), stage)
    final_column = _apply_reduction_factors(factors, column, R)
    provisional = ECPColumnReductionCertificate(column, R, stages, factors, final_column, nothing)
    verification = _ecp_column_reduction_replay_summary(provisional)
    verification.overall_ok || error("internal ECP stage certificate verification failed")
    certificate = ECPColumnReductionCertificate(column, R, stages, factors, final_column, verification)
    verify_ecp_column_reduction(certificate) || error("internal ECP stage certificate storage verification failed")
    return certificate
end

function verify_ecp_column_reduction(certificate)::Bool
    try
        replay = _ecp_column_reduction_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _ecp_column_reduction_replay_summary(certificate)
    replay = _ecp_replay_stages(certificate)
    replayed_factors = collect(replay.replayed_factors)
    replayed_final_column = _apply_reduction_factors(replayed_factors, certificate.original_column, certificate.ring)
    target = _target_reduced_column(certificate.ring, length(certificate.original_column))
    factors_match_ok = _ecp_factor_sequences_equal(certificate.factors, replayed_factors)
    final_column_ok = certificate.final_column == replayed_final_column
    target_ok = replayed_final_column == target
    overall_ok = replay.ok && factors_match_ok && final_column_ok && target_ok
    return (
        overall_ok = overall_ok,
        replay = replay,
        factors_match_ok = factors_match_ok,
        final_column_ok = final_column_ok,
        target_ok = target_ok,
        replayed_factors = replayed_factors,
        replayed_final_column = replayed_final_column,
    )
end

function _ecp_replay_stages(certificate)
    stages = certificate.stages
    length(stages) == 2 || return (; ok = false, stage_replays = (), replayed_factors = Any[])

    current_input = collect(certificate.original_column)
    stage_replays = Any[]
    replayed_factors = Any[]
    for stage in stages
        replay = _ecp_replay_stage(stage, current_input, certificate.ring)
        push!(stage_replays, replay)
        replay.ok || return (; ok = false, stage_replays = tuple(stage_replays...), replayed_factors = replayed_factors)
        append!(replayed_factors, replay.factors)
        current_input = collect(_ecp_matrix_column_to_tuple(replay.output_column))
    end

    return (; ok = true, stage_replays = tuple(stage_replays...), replayed_factors = replayed_factors)
end

function _ecp_replay_stage(stage, input_column, R)
    if stage.kind == :validation
        ok = _ecp_stage_keys_ok(stage, (:kind, :input_length, :is_unimodular)) &&
            stage.input_length == length(input_column) &&
            stage.is_unimodular === is_unimodular_column(collect(input_column), R)
        return (; ok, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
    elseif stage.kind == :unit_entry
        expected_factors = _reduce_unit_witness_column(input_column, stage.pivot_index, R)
        expected_output = _apply_reduction_factors(expected_factors, input_column, R)
        ok = _ecp_stage_keys_ok(stage, (:kind, :input_column, :pivot_index, :pivot_value, :pivot_inverse, :factors, :output_column)) &&
            stage.input_column == _ecp_column_tuple(input_column) &&
            stage.pivot_value == input_column[stage.pivot_index] &&
            stage.pivot_inverse == inv(input_column[stage.pivot_index]) &&
            _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
            stage.output_column == expected_output
        return (; ok, factors = expected_factors, output_column = expected_output)
    elseif stage.kind == :witness_unit
        witness = collect(stage.witness)
        unit_creation_factors = _witness_unit_creation_factors(input_column, witness, stage.pivot_index, R)
        created_column = collect(_ecp_matrix_column_to_tuple(_apply_reduction_factors(unit_creation_factors, input_column, R)))
        unit_replay = _ecp_replay_stage(stage.unit_stage, created_column, R)
        expected_factors = vcat(unit_replay.factors, unit_creation_factors)
        expected_output = _apply_reduction_factors(expected_factors, input_column, R)
        witness_total = sum(witness[idx] * input_column[idx] for idx in eachindex(input_column))
        ok = _ecp_stage_keys_ok(
                stage,
                (:kind, :input_column, :witness, :pivot_index, :witness_unit, :witness_unit_inverse, :unit_creation_factors, :created_column, :unit_stage, :factors, :output_column),
            ) &&
            stage.input_column == _ecp_column_tuple(input_column) &&
            length(stage.witness) == length(input_column) &&
            witness_total == one(R) &&
            stage.witness_unit == witness[stage.pivot_index] &&
            stage.witness_unit_inverse == inv(witness[stage.pivot_index]) &&
            _ecp_factor_sequences_equal(stage.unit_creation_factors, unit_creation_factors) &&
            stage.created_column == tuple(created_column...) &&
            unit_replay.ok &&
            _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
            stage.output_column == expected_output
        return (; ok, factors = expected_factors, output_column = expected_output)
    elseif stage.kind == :monicity_normalization
        ring_gens = collect(gens(R))
        variable_order = collect(stage.variable_order)
        target_variable = stage.target_variable
        target_variable_index = stage.target_variable_index
        source_variable_index = stage.source_variable_index
        source_variable = stage.source_variable
        forward_values = copy(ring_gens)
        forward_values[source_variable_index] = source_variable + stage.shift_sign * target_variable^stage.shift_power
        substituted_column = [
            _coerce_into_ring(R, evaluate(entry, forward_values), "substituted column entry")
            for entry in input_column
        ]
        transformed_replay = _ecp_replay_stage(stage.transformed_stage, substituted_column, R)
        inverse_values = copy(ring_gens)
        inverse_values[source_variable_index] = source_variable - stage.shift_sign * target_variable^stage.shift_power
        forward_substitution = _ecp_substitution_map_tuple(ring_gens, forward_values)
        inverse_substitution = _ecp_substitution_map_tuple(ring_gens, inverse_values)
        selected_index_ok = stage.selected_monic_index isa Integer &&
            1 <= stage.selected_monic_index <= length(substituted_column)
        selected_monic_entry = selected_index_ok ? substituted_column[stage.selected_monic_index] : zero(R)
        selected_monic_ok = selected_index_ok &&
            stage.selected_monic_entry == selected_monic_entry &&
            _is_monic_in_variable(selected_monic_entry, R, target_variable_index)
        first_coordinate_strategy = selected_index_ok ?
            _ecp_first_coordinate_strategy(stage.selected_monic_index) :
            :invalid
        transformed_output = _apply_reduction_factors(transformed_replay.factors, substituted_column, R)
        target = _target_reduced_column(R, length(input_column))
        inverse_substituted_factors = [
            _substitute_matrix_entries(factor, inverse_values, R)
            for factor in transformed_replay.factors
        ]
        expected_output = _apply_reduction_factors(inverse_substituted_factors, input_column, R)
        first_coordinate_move_factors = typeof(identity_matrix(R, length(input_column)))[]
        variable_change_verification = (;
            selected_monic_ok,
            transformed_reduction_ok = transformed_output == target,
            original_reduction_ok = expected_output == target,
        )
        normalized_variable_order = _ecp_normalize_variable_order(R, variable_order)
        source_order_index = findfirst(variable -> variable == source_variable, variable_order)
        target_order_index = findfirst(variable -> variable == target_variable, variable_order)
        ok = _ecp_stage_keys_ok(
                stage,
                (:kind, :input_column, :variable_order, :variable_index, :source_variable_index, :source_variable, :last_variable_index, :target_variable_index, :target_variable, :shift_power, :shift_sign, :shift_polynomial, :forward_values, :inverse_values, :forward_substitution, :inverse_substitution, :transformed_column, :selected_monic_index, :selected_monic_entry, :first_coordinate_strategy, :first_coordinate_move_factors, :first_coordinate_column, :transformed_stage, :transformed_factors, :inverse_substituted_factors, :factors, :output_column, :variable_change_verification),
            ) &&
            stage.input_column == _ecp_column_tuple(input_column) &&
            stage.variable_order == tuple(normalized_variable_order...) &&
            stage.source_variable_index == stage.variable_index &&
            stage.source_variable == ring_gens[stage.source_variable_index] &&
            stage.last_variable_index == stage.target_variable_index &&
            stage.target_variable == variable_order[end] &&
            source_order_index !== nothing &&
            target_order_index !== nothing &&
            source_order_index < target_order_index &&
            target_order_index == length(variable_order) &&
            stage.target_variable == ring_gens[stage.target_variable_index] &&
            stage.shift_polynomial == stage.shift_sign * stage.target_variable^stage.shift_power &&
            stage.forward_values == tuple(forward_values...) &&
            stage.inverse_values == tuple(inverse_values...) &&
            stage.forward_substitution == forward_substitution &&
            stage.inverse_substitution == inverse_substitution &&
            stage.transformed_column == tuple(substituted_column...) &&
            selected_monic_ok &&
            stage.first_coordinate_strategy == first_coordinate_strategy &&
            stage.first_coordinate_move_factors == first_coordinate_move_factors &&
            stage.first_coordinate_column == tuple(substituted_column...) &&
            transformed_replay.ok &&
            _ecp_factor_sequences_equal(stage.transformed_factors, transformed_replay.factors) &&
            _ecp_factor_sequences_equal(stage.inverse_substituted_factors, inverse_substituted_factors) &&
            _ecp_factor_sequences_equal(stage.factors, inverse_substituted_factors) &&
            stage.variable_change_verification == variable_change_verification &&
            variable_change_verification.selected_monic_ok &&
            variable_change_verification.transformed_reduction_ok &&
            variable_change_verification.original_reduction_ok &&
            stage.output_column == expected_output
        return (; ok, factors = inverse_substituted_factors, output_column = expected_output)
    elseif stage.kind == :embedded_three_block
        indices = collect(stage.indices)
        subcolumn = [input_column[idx] for idx in indices]
        subcertificate_ok = verify_ecp_column_reduction(stage.subcertificate) &&
            stage.subcertificate.original_column == subcolumn &&
            stage.subcertificate.ring == R
        embedded_factors = [block_embedding(factor, length(input_column), indices) for factor in stage.subcertificate.factors]
        post_block = _apply_reduction_factors(embedded_factors, input_column, R)
        pivot_idx = indices[end]
        elimination_factors = typeof(identity_matrix(R, length(input_column)))[]
        for row in 1:length(input_column)
            row == pivot_idx && continue
            coeff = -post_block[row, 1]
            coeff == zero(R) && continue
            push!(elimination_factors, elementary_matrix(length(input_column), row, pivot_idx, coeff, R))
        end
        move_factors = typeof(identity_matrix(R, length(input_column)))[]
        if pivot_idx != length(input_column)
            push!(move_factors, elementary_matrix(length(input_column), pivot_idx, length(input_column), -one(R), R))
            push!(move_factors, elementary_matrix(length(input_column), length(input_column), pivot_idx, one(R), R))
        end
        expected_factors = vcat(move_factors, elimination_factors, embedded_factors)
        expected_output = _apply_reduction_factors(expected_factors, input_column, R)
        ok = _ecp_stage_keys_ok(
                stage,
                (:kind, :input_column, :indices, :subcolumn, :subcertificate, :embedded_factors, :post_block_column, :elimination_factors, :move_factors, :factors, :output_column),
            ) &&
            stage.input_column == _ecp_column_tuple(input_column) &&
            stage.subcolumn == tuple(subcolumn...) &&
            subcertificate_ok &&
            _ecp_factor_sequences_equal(stage.embedded_factors, embedded_factors) &&
            stage.post_block_column == _ecp_matrix_column_to_tuple(post_block) &&
            _ecp_factor_sequences_equal(stage.elimination_factors, elimination_factors) &&
            _ecp_factor_sequences_equal(stage.move_factors, move_factors) &&
            _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
            stage.output_column == expected_output
        return (; ok, factors = expected_factors, output_column = expected_output)
    elseif stage.kind == :laurent_normalization
        normalization = normalize_laurent_object(input_column)
        polynomial_certificate_ok = verify_ecp_column_reduction(stage.polynomial_certificate) &&
            stage.polynomial_certificate.original_column == collect(normalization.normalized_object) &&
            stage.polynomial_certificate.ring == normalization.metadata.polynomial_ring
        lifted_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in stage.polynomial_certificate.factors]
        shift = only(normalization.metadata.shift_monomials)
        inverse_shift = only(normalization.metadata.inverse_shift_monomials)
        normalization_factors = _unit_normalization_factors(length(input_column), inverse_shift, shift, R)
        expected_factors = vcat(normalization_factors, lifted_factors)
        expected_output = _apply_reduction_factors(expected_factors, input_column, R)
        ok = _ecp_stage_keys_ok(
                stage,
                (:kind, :input_column, :normalization, :normalized_column, :polynomial_certificate, :lifted_factors, :shift, :inverse_shift, :normalization_factors, :factors, :output_column),
            ) &&
            stage.input_column == _ecp_column_tuple(input_column) &&
            stage.normalization == normalization &&
            verify_laurent_normalization(input_column, stage.normalization) &&
            stage.normalized_column == _ecp_column_tuple(normalization.normalized_object) &&
            polynomial_certificate_ok &&
            _ecp_factor_sequences_equal(stage.lifted_factors, lifted_factors) &&
            stage.shift == shift &&
            stage.inverse_shift == inverse_shift &&
            _ecp_factor_sequences_equal(stage.normalization_factors, normalization_factors) &&
            _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
            stage.output_column == expected_output
        return (; ok, factors = expected_factors, output_column = expected_output)
    end

    return (; ok = false, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
end

function _ecp_factor_sequences_equal(left, right)
    length(left) == length(right) || return false
    for idx in eachindex(left)
        left[idx] == right[idx] || return false
    end
    return true
end

function _ecp_stage_keys_ok(stage, expected)
    return propertynames(stage) == expected
end

function _ecp_column_tuple(column)
    return tuple((entry for entry in column)...)
end

function _ecp_substitution_map_tuple(variables, values)
    return tuple(((; variable = variables[idx], value = values[idx]) for idx in eachindex(variables))...)
end

function _ecp_first_monic_entry_index(column, R)
    return _ecp_first_monic_entry_index(column, R, ngens(R))
end

function _ecp_first_monic_entry_index(column, R, variable_index::Int)
    return findfirst(entry -> _is_monic_in_variable(entry, R, variable_index), column)
end

function _ecp_first_coordinate_strategy(selected_monic_index::Int)
    return selected_monic_index == 1 ? :already_first : :not_moved
end

function _ecp_matrix_column_to_tuple(column_matrix)
    ncols(column_matrix) == 1 || throw(ArgumentError("expected a column matrix"))
    return tuple((column_matrix[row, 1] for row in 1:nrows(column_matrix))...)
end

function _lift_polynomial_reduction_factor(factor, R)
    entries = [
        _coerce_into_ring(R, factor[row, col], "lifted reduction factor entry")
        for col in 1:ncols(factor), row in 1:nrows(factor)
    ]
    return matrix(R, nrows(factor), ncols(factor), vec(entries))
end

function _is_monic_in_variable(p, R, variable_index::Int)
    iszero(p) && return false
    return _leading_coefficient_in_variable(p, R, variable_index) == one(R)
end

function _is_monic_in_last_variable(p, R)
    return _is_monic_in_variable(p, R, ngens(R))
end

function _leading_coefficient_in_variable(p, R, variable_index::Int)
    variable_index < 1 && return zero(R)
    variable_index > ngens(R) && return zero(R)

    target_degree = degree(p, variable_index)
    target_degree < 0 && return zero(R)

    ring_gens = collect(gens(R))
    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        exponents[variable_index] == target_degree || continue
        term = R(coeff)
        for idx in eachindex(ring_gens)
            idx == variable_index && continue
            exponent = exponents[idx]
            exponent == 0 && continue
            term *= ring_gens[idx]^exponent
        end
        total += term
    end

    return total
end

function _leading_coefficient_in_last_variable(p, R)
    return _leading_coefficient_in_variable(p, R, ngens(R))
end

function _substitute_matrix_entries(M, values::AbstractVector, R)
    entries = [_coerce_into_ring(R, evaluate(M[row, col], values), "substituted matrix entry") for col in 1:ncols(M), row in 1:nrows(M)]
    return matrix(R, nrows(M), ncols(M), vec(entries))
end

function _unit_normalization_factors(n::Int, u, uinv, R)
    factors = typeof(identity_matrix(R, n))[]
    u == one(R) && return factors

    i = n - 1
    j = n

    for (row, col, coeff) in (
        (i, j, u - one(R)),
        (j, i, one(R)),
        (i, j, uinv - one(R)),
        (j, i, -u),
    )
        push!(factors, elementary_matrix(n, row, col, coeff, R))
    end

    return factors
end

function _throw_unsupported_unimodular_column_reduction(column::AbstractVector, R)
    n = length(column)
    profile = _is_laurent_polynomial_ring(R) ? "Laurent-normalized" : "ordinary polynomial"
    throw(ArgumentError("unsupported exact unimodular column reduction for $(profile) column of length $(n): no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies"))
end
