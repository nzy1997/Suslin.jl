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

struct ECPLinkWitnessRecord
    original_column
    ring
    variable_order
    selected_variable_index::Int
    selected_variable
    selected_monic_index::Int
    selected_monic_entry
    residue_probes
    tail_reductions
    resultants
    bezout_coefficients
    coverage_multipliers
    path_points
    metadata
    verification
end

struct ECPLinkStepCertificate
    original_column
    ring
    link_witness::ECPLinkWitnessRecord
    path_points
    path_columns
    segments
    lower_variable_column
    transformed_column
    forward_factors::Vector
    reduction_factors::Vector
    verification
end

struct ECPInductionNormalityCertificate
    original_column
    ring
    link_step::ECPLinkStepCertificate
    lower_variable_column
    lower_reduction_certificate
    lower_variable_factors::Vector
    lifted_lower_variable_factors::Vector
    normality_witness
    normality_rewrite
    final_factors::Vector
    final_column
    verification
end

struct ECPStagedColumnReductionCertificate
    original_column
    ring
    monicity
    link_step::ECPLinkStepCertificate
    lower_reduction
    normality_witness
    induction_normality::ECPInductionNormalityCertificate
    factors::Vector
    final_column
    verification
end

function reduce_unimodular_column(v::AbstractVector, R)
    column = _validated_unimodular_column(v, R)
    staged = _ecp_public_staged_reduction_certificate(column, R)
    staged !== nothing && return staged.factors
    return _ecp_column_reduction_certificate_validated(column, R).factors
end

function ecp_column_reduction_certificate(v::AbstractVector, R)
    column = _validated_unimodular_column(v, R)
    return _ecp_column_reduction_certificate_validated(column, R)
end

function diagnose_unimodular_column_reduction(v::AbstractVector, R)
    ring_profile = _column_reduction_ring_profile(R)
    column_length = length(v)
    validation = _diagnose_unimodular_column_preconditions(v, R, ring_profile, column_length)
    validation.status == :ok || return validation.diagnostic
    return _diagnose_unimodular_column_reduction_validated(
        validation.column,
        R,
        ring_profile,
        column_length,
    )
end

function _column_reduction_diagnostic(
    status::Symbol,
    failure_code,
    ring_profile,
    column_length::Int,
    attempted_stages,
    message::AbstractString,
)
    return (;
        status,
        failure_code,
        ring_profile,
        column_length,
        attempted_stages = tuple(attempted_stages...),
        message = String(message),
    )
end

function _column_reduction_ring_profile(R)
    kind = try
        _is_laurent_polynomial_ring(R) ? :laurent_polynomial : :polynomial
    catch err
        err isa InterruptException && rethrow()
        :unknown
    end
    coefficient_ring = try
        string(base_ring(R))
    catch err
        err isa InterruptException && rethrow()
        ""
    end
    generator_names = try
        tuple(string.(gens(R))...)
    catch err
        err isa InterruptException && rethrow()
        ()
    end
    return (; kind, coefficient_ring, generators = generator_names)
end

function _diagnose_unimodular_column_preconditions(v::AbstractVector, R, ring_profile, column_length::Int)
    try
        Base.require_one_based_indexing(v)
    catch err
        err isa InterruptException && rethrow()
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :unsupported_indexing,
                ring_profile,
                column_length,
                Symbol[],
                "v must use one-based indexing",
            ),
        )
    end

    if column_length < 3
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :column_too_short,
                ring_profile,
                column_length,
                Symbol[],
                "v must have length at least 3",
            ),
        )
    end

    column = try
        [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:column_length]
    catch err
        err isa InterruptException && rethrow()
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :column_not_over_ring,
                ring_profile,
                column_length,
                Symbol[],
                _column_reduction_error_message(err),
            ),
        )
    end

    is_unimodular = try
        is_unimodular_column(column, R)
    catch err
        err isa InterruptException && rethrow()
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :unimodularity_check_failed,
                ring_profile,
                column_length,
                Symbol[],
                _column_reduction_error_message(err),
            ),
        )
    end

    is_unimodular || return (;
        status = :precondition_failed,
        diagnostic = _column_reduction_diagnostic(
            :precondition_failed,
            :not_unimodular,
            ring_profile,
            column_length,
            Symbol[],
            "v must be a unimodular column",
        ),
    )
    return (; status = :ok, column)
end

function _diagnose_unimodular_column_reduction_validated(column::AbstractVector, R, ring_profile, column_length::Int)
    attempted = Symbol[]
    result = _is_laurent_polynomial_ring(R) ?
        _diagnose_laurent_unimodular_column_reduction(column, R, attempted) :
        _diagnose_polynomial_unimodular_column_reduction(column, R, attempted)

    if result.supported
        return _column_reduction_diagnostic(
            :supported,
            nothing,
            ring_profile,
            column_length,
            attempted,
            "exact unimodular column reduction is supported by $(result.stage)",
        )
    end

    failure_code = _is_laurent_polynomial_ring(R) ?
        :unsupported_laurent_column_family :
        :unsupported_polynomial_column_family
    return _column_reduction_diagnostic(
        :unsupported,
        failure_code,
        ring_profile,
        column_length,
        attempted,
        _unsupported_unimodular_column_reduction_message(column, R),
    )
end

function _ecp_column_reduction_certificate_validated(column::AbstractVector, R)
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

function ecp_link_witness(
    v::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    selected_monic_index::Integer = 1,
    supplied_link_witness = nothing,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP link witnesses currently support ordinary polynomial columns only"))
    supplied_link_witness === nothing &&
        throw(ArgumentError("Park-Woodburn ECP link witness extraction is not implemented; pass supplied_link_witness metadata with source = :supplied_link_witness"))

    column = _validated_unimodular_column(v, R)
    normalized_order = _ecp_normalize_variable_order(R, variable_order)
    isempty(normalized_order) && selected_variable === nothing &&
        throw(ArgumentError("variable_order must contain at least one generator when selected_variable is not provided"))
    selected_variable = selected_variable === nothing ? first(normalized_order) : selected_variable
    selected_variable_index = _ecp_selected_variable_index(R, selected_variable)
    selected_monic_index = Int(selected_monic_index)
    selected_monic_index == 1 ||
        throw(ArgumentError("Park-Woodburn ECP link witnesses require the selected monic entry to be first"))
    _is_monic_in_variable(column[selected_monic_index], R, selected_variable_index) ||
        throw(ArgumentError("selected first entry must be monic in the selected variable"))

    metadata = (; source = _ecp_link_field(supplied_link_witness, :source))
    metadata.source == :supplied_link_witness ||
        throw(ArgumentError("supplied ECP link witness metadata must use source = :supplied_link_witness"))
    record, verification = try
        local replay_record = ECPLinkWitnessRecord(
            tuple(column...),
            R,
            tuple(normalized_order...),
            selected_variable_index,
            gens(R)[selected_variable_index],
            selected_monic_index,
            column[selected_monic_index],
            tuple(_ecp_link_field(supplied_link_witness, :residue_probes)...),
            tuple(_ecp_link_field(supplied_link_witness, :tail_reductions)...),
            tuple(_ecp_link_field(supplied_link_witness, :resultants)...),
            tuple(_ecp_link_field(supplied_link_witness, :bezout_coefficients)...),
            tuple(_ecp_link_field(supplied_link_witness, :coverage_multipliers)...),
            tuple(_ecp_link_field(supplied_link_witness, :path_points)...),
            metadata,
            nothing,
        )
        replay_record, _ecp_link_witness_replay_summary(replay_record)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError("supplied Park-Woodburn ECP link witness data failed exact replay verification"))
    end
    verification.overall_ok ||
        throw(ArgumentError("supplied Park-Woodburn ECP link witness data failed exact replay verification"))
    stored = ECPLinkWitnessRecord(
        record.original_column,
        record.ring,
        record.variable_order,
        record.selected_variable_index,
        record.selected_variable,
        record.selected_monic_index,
        record.selected_monic_entry,
        record.residue_probes,
        record.tail_reductions,
        record.resultants,
        record.bezout_coefficients,
        record.coverage_multipliers,
        record.path_points,
        record.metadata,
        verification,
    )
    verify_ecp_link_witness(stored) ||
        throw(ArgumentError("stored Park-Woodburn ECP link witness data failed exact replay verification"))
    return stored
end

function ecp_link_step_certificate(
    v::AbstractVector,
    R;
    link_witness = nothing,
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    selected_monic_index::Integer = 1,
    supplied_link_witness = nothing,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP link steps currently support ordinary polynomial columns only"))
    witness = link_witness === nothing ?
        ecp_link_witness(v, R; variable_order, selected_variable, selected_monic_index, supplied_link_witness) :
        link_witness
    verify_ecp_link_witness(witness) ||
        throw(ArgumentError("ECP link step requires a verified Park-Woodburn link witness"))
    _same_base_ring(witness.ring, R) ||
        throw(ArgumentError("ECP link step input ring must match the link witness ring"))

    column = _validated_unimodular_column(v, R)
    tuple(column...) == witness.original_column ||
        throw(ArgumentError("ECP link step input column must match the link witness column"))

    path_columns = _ecp_link_step_path_columns(witness)
    segments = _ecp_link_step_segments(witness, path_columns)
    forward_factors = _ecp_link_step_forward_factors(segments)
    reduction_factors = _ecp_link_step_reduction_factors(segments)
    provisional = ECPLinkStepCertificate(
        tuple(column...),
        R,
        witness,
        witness.path_points,
        path_columns,
        segments,
        first(path_columns),
        tuple(column...),
        forward_factors,
        reduction_factors,
        nothing,
    )
    verification = _ecp_link_step_replay_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("constructed Park-Woodburn ECP link step failed exact replay verification"))
    certificate = ECPLinkStepCertificate(
        provisional.original_column,
        provisional.ring,
        provisional.link_witness,
        provisional.path_points,
        provisional.path_columns,
        provisional.segments,
        provisional.lower_variable_column,
        provisional.transformed_column,
        provisional.forward_factors,
        provisional.reduction_factors,
        verification,
    )
    verify_ecp_link_step_certificate(certificate) ||
        throw(ArgumentError("stored Park-Woodburn ECP link step failed exact replay verification"))
    return certificate
end

function ecp_induction_normality_certificate(
    v::AbstractVector,
    R;
    link_step = nothing,
    lower_reduction = nothing,
    normality_witness = nothing,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP induction/normality currently supports ordinary polynomial columns only"))
    link_step === nothing &&
        throw(ArgumentError("ECP induction/normality requires a verified link-step certificate"))
    verify_ecp_link_step_certificate(link_step) ||
        throw(ArgumentError("ECP induction/normality requires a verified link-step certificate"))
    _same_base_ring(link_step.ring, R) ||
        throw(ArgumentError("ECP induction/normality input ring must match the link-step ring"))

    column = _validated_unimodular_column(v, R)
    tuple(column...) == link_step.original_column ||
        throw(ArgumentError("ECP induction/normality input column must match the link-step column"))

    lower_column = collect(link_step.lower_variable_column)
    lower_certificate, lower_factors = _ecp_verified_lower_reduction(lower_reduction, lower_column, R)
    lifted_lower_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in lower_factors]
    normality_rewrite = _ecp_induction_normality_rewrite(
        normality_witness,
        lower_column,
        lifted_lower_factors,
        R,
    )
    final_factors = vcat(lifted_lower_factors, normality_rewrite.rewrite_factors, link_step.reduction_factors)
    final_column = _apply_reduction_factors(final_factors, column, R)
    provisional = ECPInductionNormalityCertificate(
        tuple(column...),
        R,
        link_step,
        tuple(lower_column...),
        lower_certificate,
        lower_factors,
        lifted_lower_factors,
        normality_witness,
        normality_rewrite,
        final_factors,
        final_column,
        nothing,
    )
    verification = _ecp_induction_normality_replay_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("constructed ECP induction/normality certificate failed exact replay verification"))
    certificate = ECPInductionNormalityCertificate(
        provisional.original_column,
        provisional.ring,
        provisional.link_step,
        provisional.lower_variable_column,
        provisional.lower_reduction_certificate,
        provisional.lower_variable_factors,
        provisional.lifted_lower_variable_factors,
        provisional.normality_witness,
        provisional.normality_rewrite,
        provisional.final_factors,
        provisional.final_column,
        verification,
    )
    verify_ecp_induction_normality_certificate(certificate) ||
        throw(ArgumentError("stored ECP induction/normality certificate failed exact replay verification"))
    return certificate
end

function ecp_staged_column_reduction_certificate(
    v::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    selected_monic_index::Integer = 1,
    supplied_link_witness = nothing,
    lower_reduction = nothing,
    normality_witness = nothing,
)
    column = _validated_unimodular_column(v, R)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP staged public column pipeline currently supports ordinary polynomial columns only"))
    normalized_order = _ecp_normalize_variable_order(R, variable_order)
    isempty(normalized_order) && selected_variable === nothing &&
        throw(ArgumentError("variable_order must contain at least one generator when selected_variable is not provided"))
    selected_variable = selected_variable === nothing ? first(normalized_order) : selected_variable
    link_witness = supplied_link_witness === nothing ?
        _ecp_default_public_link_witness(column, R, selected_variable) :
        supplied_link_witness
    link = ecp_link_step_certificate(
        column,
        R;
        variable_order,
        selected_variable,
        selected_monic_index,
        supplied_link_witness = link_witness,
    )
    lower_column = collect(link.lower_variable_column)
    lower_certificate, lower_factors = _ecp_verified_lower_reduction(lower_reduction, lower_column, R)
    lower_reduction_input = lower_certificate === nothing ? lower_factors : lower_certificate
    normality = normality_witness === nothing ?
        _ecp_default_public_normality_witness(lower_factors, length(lower_column), R) :
        normality_witness
    induction = ecp_induction_normality_certificate(
        column,
        R;
        link_step = link,
        lower_reduction = lower_reduction_input,
        normality_witness = normality,
    )
    monicity = (;
        source = :link_witness,
        variable_order = link.link_witness.variable_order,
        selected_variable_index = link.link_witness.selected_variable_index,
        selected_variable = link.link_witness.selected_variable,
        selected_monic_index = link.link_witness.selected_monic_index,
        selected_monic_entry = link.link_witness.selected_monic_entry,
        selected_monic_ok = link.link_witness.verification.selected_monic_ok,
    )
    factors = collect(induction.final_factors)
    final_column = _apply_reduction_factors(factors, column, R)
    provisional = ECPStagedColumnReductionCertificate(
        tuple(column...),
        R,
        monicity,
        link,
        lower_reduction_input,
        normality,
        induction,
        factors,
        final_column,
        nothing,
    )
    verification = _ecp_staged_column_reduction_replay_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("constructed ECP staged public column pipeline failed exact replay verification"))
    certificate = ECPStagedColumnReductionCertificate(
        provisional.original_column,
        provisional.ring,
        provisional.monicity,
        provisional.link_step,
        provisional.lower_reduction,
        provisional.normality_witness,
        provisional.induction_normality,
        provisional.factors,
        provisional.final_column,
        verification,
    )
    verify_ecp_staged_column_reduction(certificate) ||
        throw(ArgumentError("constructed ECP staged public column pipeline failed exact replay verification"))
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

function _laurent_witness_solve_failure(err)::Bool
    err isa ErrorException || return false
    message = sprint(showerror, err)
    return occursin("No exact solution exists for A * U = B", message) ||
        occursin("not liftable to the given generating system", message)
end

function _laurent_unimodular_witness(column::AbstractVector, R; solver = solve_laurent_linear)
    _is_laurent_polynomial_ring(R) || return nothing
    row = matrix(R, 1, length(column), collect(column))
    rhs = matrix(R, 1, 1, [one(R)])

    solution = try
        solver(row, rhs)
    catch err
        err isa InterruptException && rethrow()
        _laurent_witness_solve_failure(err) && return nothing
        rethrow()
    end

    row * solution == rhs || return nothing
    return [solution[idx, 1] for idx in 1:nrows(solution)]
end

function _reduce_via_laurent_witness_unit_certificate(column::AbstractVector, R)
    witness = _laurent_unimodular_witness(column, R)
    witness === nothing && return nothing

    witness_unit_idx = findfirst(is_unit, witness)
    witness_unit_idx === nothing && return nothing

    return _witness_unit_reduction_certificate_stage(column, witness, witness_unit_idx, R)
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

function _ecp_link_field(source, field::Symbol)
    hasproperty(source, field) || throw(ArgumentError("supplied ECP link witness missing field $(field)"))
    return getproperty(source, field)
end

function _ecp_selected_variable_index(R, variable)
    ring_gens = collect(gens(R))
    idx = _ecp_variable_order_match_index(ring_gens, variable)
    idx === nothing && throw(ArgumentError("selected_variable must be a generator of R"))
    return idx
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

    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return unit_creation

    witness_unit = _reduce_via_laurent_witness_unit_certificate(column, R)
    witness_unit !== nothing && return witness_unit

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

function verify_ecp_link_witness(record)::Bool
    try
        replay = _ecp_link_witness_replay_summary(record)
        return replay.overall_ok && record.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function verify_ecp_link_step_certificate(certificate)::Bool
    try
        replay = _ecp_link_step_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function verify_ecp_induction_normality_certificate(certificate)::Bool
    try
        replay = _ecp_induction_normality_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function verify_ecp_staged_column_reduction(certificate)::Bool
    try
        replay = _ecp_staged_column_reduction_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _ecp_verified_lower_reduction(lower_reduction, lower_column, R)
    if lower_reduction === nothing
        certificate = ecp_column_reduction_certificate(lower_column, R)
        return certificate, certificate.factors
    end

    if lower_reduction isa ECPColumnReductionCertificate
        verify_ecp_column_reduction(lower_reduction) ||
            throw(ArgumentError("lower-variable reduction certificate does not verify"))
        _same_base_ring(lower_reduction.ring, R) ||
            throw(ArgumentError("lower-variable reduction ring must match the original ring"))
        lower_reduction.original_column == lower_column ||
            throw(ArgumentError("lower-variable reduction column must match v(0)"))
        return lower_reduction, lower_reduction.factors
    end

    _ecp_is_concrete_factor_sequence(lower_reduction) ||
        throw(ArgumentError("lower-variable reduction must be an ECP certificate or a concrete elementary factor sequence"))
    factors = collect(lower_reduction)
    n = length(lower_column)
    all(factor -> _ecp_is_elementary_factor(factor, R, n), factors) ||
        throw(ArgumentError("lower-variable factor sequence must contain only elementary factors over R"))
    _apply_reduction_factors(factors, lower_column, R) == _target_reduced_column(R, length(lower_column)) ||
        throw(ArgumentError("lower-variable factor sequence does not reduce v(0) to e_n"))
    return nothing, factors
end

function _ecp_is_concrete_factor_sequence(value)
    value isa AbstractVector && return true
    value isa Tuple && !(value isa NamedTuple) && return true
    return false
end

function _ecp_induction_normality_rewrite(normality_witness, lower_column, lifted_lower_factors, R)
    normality_witness === nothing &&
        throw(ArgumentError("ECP induction/normality requires explicit normality witness data"))
    _ecp_normality_witness_keys_ok(normality_witness) ||
        throw(ArgumentError("normality witness must contain source, conjugator, sl2_indices, and sl2_entry"))
    normality_witness.source == :supplied_normality_witness ||
        throw(ArgumentError("normality witness must use source = :supplied_normality_witness"))

    n = length(lower_column)
    lower_product = _factor_sequence_product(lifted_lower_factors, R, n)
    conjugator = normality_witness.conjugator
    nrows(conjugator) == n && ncols(conjugator) == n && _same_base_ring(base_ring(conjugator), R) ||
        throw(ArgumentError("normality witness conjugator must be an n by n matrix over R"))
    conjugator * lower_product == identity_matrix(R, n) ||
        throw(ArgumentError("normality witness conjugator must invert the lifted lower-variable reduction"))

    sl2_indices = tuple(normality_witness.sl2_indices...)
    length(sl2_indices) == 2 || throw(ArgumentError("normality witness sl2_indices must contain two indices"))
    fixed_index = Int(sl2_indices[1])
    moving_index = Int(sl2_indices[2])
    fixed_index == n || throw(ArgumentError("normality witness must use e_n as the fixed SL_2 coordinate"))
    1 <= moving_index <= n && moving_index != fixed_index ||
        throw(ArgumentError("normality witness moving index must be distinct and in range"))
    entry = _coerce_into_ring(R, normality_witness.sl2_entry, "normality witness sl2_entry")
    entry != zero(R) || throw(ArgumentError("normality witness must record a non-identity SL_2 contribution"))

    sl2_block = elementary_matrix(2, 1, 2, entry, R)
    sl2_embedding = block_embedding(sl2_block, n, sl2_indices)
    rewrite_factors = try
        realize_conjugate_elementary(conjugator, fixed_index, moving_index, entry)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("normality witness could not be rewritten into elementary factors"))
    end
    rewrite_product = _factor_sequence_product(rewrite_factors, R, n)
    expected_rewrite_product = conjugator * sl2_embedding * lower_product
    lower_matrix = matrix(R, n, 1, collect(lower_column))
    fixed_lower_column_ok = rewrite_product * lower_matrix == lower_matrix
    rewrite_product_ok = rewrite_product == expected_rewrite_product
    return (;
        source = :supplied_normality_witness,
        conjugator,
        lower_product,
        sl2_indices,
        fixed_index,
        moving_index,
        sl2_entry = entry,
        sl2_block,
        sl2_embedding,
        rewrite_factors,
        rewrite_product,
        expected_rewrite_product,
        rewrite_product_ok,
        fixed_lower_column_ok,
        overall_ok = rewrite_product_ok && fixed_lower_column_ok,
    )
end

function _ecp_normality_witness_keys_ok(normality_witness)
    names = propertynames(normality_witness)
    return all(field -> field in names, (:source, :conjugator, :sl2_indices, :sl2_entry))
end

function _ecp_induction_normality_replay_summary(certificate)
    R = certificate.ring
    n = length(certificate.original_column)
    link_step_ok = verify_ecp_link_step_certificate(certificate.link_step)
    input_ok = link_step_ok && certificate.original_column == certificate.link_step.original_column
    lower_variable_column_ok = link_step_ok &&
        certificate.lower_variable_column == certificate.link_step.lower_variable_column

    lower_certificate, lower_factors = try
        _ecp_verified_lower_reduction(
            certificate.lower_reduction_certificate === nothing ?
                certificate.lower_variable_factors :
                certificate.lower_reduction_certificate,
            collect(certificate.lower_variable_column),
            R,
        )
    catch err
        err isa InterruptException && rethrow()
        nothing, Any[]
    end
    lower_reduction_ok = _ecp_factor_sequences_equal(certificate.lower_variable_factors, lower_factors)
    lifted_lower_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in lower_factors]
    lifted_lower_factors_ok = _ecp_factor_sequences_equal(certificate.lifted_lower_variable_factors, lifted_lower_factors)

    normality_rewrite = try
        _ecp_induction_normality_rewrite(
            certificate.normality_witness,
            collect(certificate.lower_variable_column),
            certificate.lifted_lower_variable_factors,
            R,
        )
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    normality_rewrite_ok = normality_rewrite !== nothing &&
        certificate.normality_rewrite == normality_rewrite &&
        normality_rewrite.overall_ok

    expected_final_factors = normality_rewrite === nothing ?
        Any[] :
        vcat(certificate.lifted_lower_variable_factors, normality_rewrite.rewrite_factors, certificate.link_step.reduction_factors)
    final_factors_ok = _ecp_factor_sequences_equal(certificate.final_factors, expected_final_factors)
    final_column = final_factors_ok ?
        _apply_reduction_factors(certificate.final_factors, collect(certificate.original_column), R) :
        zero_matrix(R, n, 1)
    final_column_ok = certificate.final_column == final_column
    final_reduction_ok = final_column == _target_reduced_column(R, n)
    final_factors_elementary_ok = all(factor -> _ecp_is_elementary_factor(factor, R, n), certificate.final_factors)
    overall_ok = link_step_ok &&
        input_ok &&
        lower_variable_column_ok &&
        lower_reduction_ok &&
        lifted_lower_factors_ok &&
        normality_rewrite_ok &&
        final_factors_ok &&
        final_column_ok &&
        final_reduction_ok &&
        final_factors_elementary_ok
    return (;
        overall_ok,
        link_step_ok,
        input_ok,
        lower_variable_column_ok,
        lower_reduction_ok,
        lifted_lower_factors_ok,
        normality_rewrite_ok,
        final_factors_ok,
        final_column_ok,
        final_reduction_ok,
        final_factors_elementary_ok,
    )
end

function _ecp_is_elementary_factor(factor, R, n::Int)
    try
        nrows(factor) == n || return false
        ncols(factor) == n || return false
        _same_base_ring(base_ring(factor), R) || return false
        identity = identity_matrix(R, n)
        positions = [(row, col) for row in 1:n, col in 1:n if factor[row, col] != identity[row, col]]
        length(positions) == 1 || return false
        row, col = only(positions)
        return row != col
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

_ecp_namedtuple_keys_exact(value, expected::Tuple) = Tuple(propertynames(value)) == expected

function _ecp_generator_belongs_to_ring(R, generator)
    return any(ring_generator -> ring_generator == generator, gens(R))
end

function _ecp_replayed_residue_probe(probe)
    return (;
        id = probe.id,
        kind = probe.kind,
        maximal_ideal_generators = tuple(probe.maximal_ideal_generators...),
    )
end

function _ecp_link_witness_replay_summary(record)
    R = record.ring
    v1 = record.original_column[1]
    tail_entries = record.original_column[2:end]
    ring_generators = gens(R)
    expected_probe_fields = (:id, :kind, :maximal_ideal_generators)

    metadata_ok = record.metadata == (; source = :supplied_link_witness)
    normalized_variable_order = try
        tuple(_ecp_normalize_variable_order(R, record.variable_order)...)
    catch
        ()
    end
    variable_order_ok = normalized_variable_order == record.variable_order
    selected_variable_index_ok = 1 <= record.selected_variable_index <= length(ring_generators)
    replayed_selected_variable = selected_variable_index_ok ? ring_generators[record.selected_variable_index] : nothing
    selected_variable_ok = selected_variable_index_ok &&
        record.selected_variable == replayed_selected_variable
    selected_variable_order_ok = variable_order_ok &&
        count(==(record.selected_variable), record.variable_order) == 1
    selected_monic_ok = selected_variable_index_ok &&
        record.selected_monic_index == 1 &&
        record.selected_monic_entry == v1 &&
        _is_monic_in_variable(record.selected_monic_entry, R, record.selected_variable_index)
    probe_shape_ok = all(record.residue_probes) do probe
        _ecp_namedtuple_keys_exact(probe, expected_probe_fields) || return false
        generators = tuple(probe.maximal_ideal_generators...)
        return all(generator -> _ecp_generator_belongs_to_ring(R, generator), generators)
    end
    replayed_residue_probes = probe_shape_ok ?
        tuple((_ecp_replayed_residue_probe(probe) for probe in record.residue_probes)...) :
        ()
    stored_probe_metadata_ok = if isnothing(record.verification)
        true
    elseif hasproperty(record.verification, :replayed_residue_probes)
        record.verification.replayed_residue_probes == replayed_residue_probes
    else
        false
    end
    metadata_shape_ok = probe_shape_ok && stored_probe_metadata_ok

    probe_ids = tuple((probe.id for probe in record.residue_probes)...)
    lengths_ok = length(record.tail_reductions) == length(record.residue_probes) &&
        length(record.resultants) == length(record.tail_reductions) &&
        length(record.bezout_coefficients) == length(record.tail_reductions) &&
        length(record.coverage_multipliers) == length(record.tail_reductions) &&
        length(record.path_points) == length(record.tail_reductions) + 1

    recomputed_tail_tilde_Gs = Any[]
    tail_reduction_ok = lengths_ok && metadata_shape_ok
    for (idx, tail) in enumerate(record.tail_reductions)
        hasproperty(tail, :probe_id) || return _ecp_invalid_link_replay(record, metadata_ok, variable_order_ok, selected_variable_index_ok, selected_variable_ok, selected_variable_order_ok, selected_monic_ok, metadata_shape_ok, lengths_ok, recomputed_tail_tilde_Gs)
        hasproperty(tail, :G) || return _ecp_invalid_link_replay(record, metadata_ok, variable_order_ok, selected_variable_index_ok, selected_variable_ok, selected_variable_order_ok, selected_monic_ok, metadata_shape_ok, lengths_ok, recomputed_tail_tilde_Gs)
        hasproperty(tail, :lifted_tail_coefficients) || return _ecp_invalid_link_replay(record, metadata_ok, variable_order_ok, selected_variable_index_ok, selected_variable_ok, selected_variable_order_ok, selected_monic_ok, metadata_shape_ok, lengths_ok, recomputed_tail_tilde_Gs)
        hasproperty(tail, :tilde_G) || return _ecp_invalid_link_replay(record, metadata_ok, variable_order_ok, selected_variable_index_ok, selected_variable_ok, selected_variable_order_ok, selected_monic_ok, metadata_shape_ok, lengths_ok, recomputed_tail_tilde_Gs)
        coeffs = tuple(tail.lifted_tail_coefficients...)
        length(coeffs) == length(tail_entries) || (tail_reduction_ok = false; push!(recomputed_tail_tilde_Gs, zero(R)); continue)
        recomputed_tilde_G = zero(R)
        for j in eachindex(coeffs)
            recomputed_tilde_G += _coerce_into_ring(R, coeffs[j], "lifted tail coefficient") * tail_entries[j]
        end
        push!(recomputed_tail_tilde_Gs, recomputed_tilde_G)
        tail_reduction_ok &= tail.probe_id in probe_ids
        tail_reduction_ok &= tail.probe_id == record.residue_probes[idx].id
        tail_reduction_ok &= tail.G == recomputed_tilde_G
        tail_reduction_ok &= tail.tilde_G == recomputed_tilde_G
    end
    recomputed_tail_tilde_Gs_tuple = tuple(recomputed_tail_tilde_Gs...)

    recomputed_resultants = Any[]
    resultants_ok = lengths_ok && tail_reduction_ok
    for (idx, tilde_G) in enumerate(recomputed_tail_tilde_Gs_tuple)
        recomputed_resultant = resultant(v1, tilde_G, record.selected_variable_index)
        push!(recomputed_resultants, recomputed_resultant)
        resultants_ok &= record.resultants[idx] == recomputed_resultant
    end
    recomputed_resultants_tuple = tuple(recomputed_resultants...)

    bezout_ok = lengths_ok && tail_reduction_ok && resultants_ok
    for idx in eachindex(record.bezout_coefficients)
        bezout = record.bezout_coefficients[idx]
        hasproperty(bezout, :f) || return _ecp_invalid_link_replay(record, metadata_ok, variable_order_ok, selected_variable_index_ok, selected_variable_ok, selected_variable_order_ok, selected_monic_ok, metadata_shape_ok, lengths_ok, recomputed_tail_tilde_Gs_tuple, recomputed_resultants_tuple)
        hasproperty(bezout, :h) || return _ecp_invalid_link_replay(record, metadata_ok, variable_order_ok, selected_variable_index_ok, selected_variable_ok, selected_variable_order_ok, selected_monic_ok, metadata_shape_ok, lengths_ok, recomputed_tail_tilde_Gs_tuple, recomputed_resultants_tuple)
        bezout_identity = bezout.f * v1 + bezout.h * recomputed_tail_tilde_Gs_tuple[idx]
        bezout_ok &= bezout_identity == recomputed_resultants_tuple[idx]
    end

    coverage_total = zero(R)
    coverage_ok = lengths_ok && resultants_ok
    for idx in eachindex(record.coverage_multipliers)
        coverage_total += record.coverage_multipliers[idx] * recomputed_resultants_tuple[idx]
    end
    coverage_ok &= coverage_total == one(R)

    path_ok = lengths_ok && coverage_ok
    if path_ok
        path_ok &= record.path_points[1] == zero(R)
        for idx in eachindex(record.coverage_multipliers)
            expected_step = record.coverage_multipliers[idx] * recomputed_resultants_tuple[idx] * record.selected_variable
            actual_step = record.path_points[idx + 1] - record.path_points[idx]
            path_ok &= actual_step == expected_step
        end
        path_ok &= record.path_points[end] == record.selected_variable
    end

    overall_ok = metadata_ok &&
        variable_order_ok &&
        selected_variable_index_ok &&
        selected_variable_ok &&
        selected_variable_order_ok &&
        selected_monic_ok &&
        metadata_shape_ok &&
        lengths_ok &&
        tail_reduction_ok &&
        resultants_ok &&
        bezout_ok &&
        coverage_ok &&
        path_ok
    return (;
        overall_ok,
        metadata_ok,
        variable_order_ok,
        normalized_variable_order,
        selected_variable_index_ok,
        selected_variable_ok,
        replayed_selected_variable,
        selected_variable_order_ok,
        selected_monic_ok,
        metadata_shape_ok,
        lengths_ok,
        replayed_residue_probes,
        tail_reduction_ok,
        resultants_ok,
        bezout_ok,
        coverage_ok,
        path_ok,
        recomputed_tail_tilde_Gs = recomputed_tail_tilde_Gs_tuple,
        recomputed_resultants = recomputed_resultants_tuple,
        coverage_total,
    )
end

function _ecp_invalid_link_replay(
    record,
    metadata_ok,
    variable_order_ok,
    selected_variable_index_ok,
    selected_variable_ok,
    selected_variable_order_ok,
    selected_monic_ok,
    metadata_shape_ok,
    lengths_ok,
    recomputed_tail_tilde_Gs = (),
    recomputed_resultants = (),
)
    return (;
        overall_ok = false,
        metadata_ok,
        variable_order_ok,
        normalized_variable_order = (),
        selected_variable_index_ok,
        selected_variable_ok,
        replayed_selected_variable = nothing,
        selected_variable_order_ok,
        selected_monic_ok,
        metadata_shape_ok,
        lengths_ok,
        replayed_residue_probes = (),
        tail_reduction_ok = false,
        resultants_ok = false,
        bezout_ok = false,
        coverage_ok = false,
        path_ok = false,
        recomputed_tail_tilde_Gs = tuple(recomputed_tail_tilde_Gs...),
        recomputed_resultants = tuple(recomputed_resultants...),
        coverage_total = zero(record.ring),
    )
end

function _ecp_link_step_path_columns(witness::ECPLinkWitnessRecord)
    return tuple((
        _ecp_evaluate_at_selected_variable(
            witness.original_column,
            witness.selected_variable_index,
            point,
            witness.ring,
        )
        for point in witness.path_points
    )...)
end

function _ecp_evaluate_at_selected_variable(values, selected_variable_index::Int, point, R)
    substitutions = collect(gens(R))
    1 <= selected_variable_index <= length(substitutions) ||
        throw(ArgumentError("selected variable index must be a generator index"))
    substitutions[selected_variable_index] = _coerce_into_ring(R, point, "path point")
    return tuple((
        _coerce_into_ring(R, evaluate(value, substitutions), "path column entry")
        for value in values
    )...)
end

function _ecp_link_step_segments(witness::ECPLinkWitnessRecord, path_columns)
    return tuple((
        _ecp_link_step_segment(witness, idx, path_columns)
        for idx in eachindex(witness.tail_reductions)
    )...)
end

function _ecp_link_step_segment(witness::ECPLinkWitnessRecord, idx::Int, path_columns)
    R = witness.ring
    n = length(witness.original_column)
    from_path_point = witness.path_points[idx]
    to_path_point = witness.path_points[idx + 1]
    delta = to_path_point - from_path_point
    from_column = path_columns[idx]
    to_column = path_columns[idx + 1]
    link_identity = _ecp_link_step_identity(witness, idx, from_column, to_column, delta)
    support_family = _ecp_link_step_supported_family(witness)
    transport = _ecp_link_step_identity_transport(R, from_column, to_column, link_identity, support_family)

    sl2_block = identity_matrix(R, 2)
    sl2_embedding = block_embedding(sl2_block, n, (1, 2))
    forward_factors = copy(transport.factors)
    inverse_factors = _ecp_inverse_factor_sequence(forward_factors)
    elementary_factors = copy(forward_factors)
    verification = _ecp_link_step_segment_verification(
        R,
        from_column,
        to_column,
        sl2_block,
        sl2_embedding,
        elementary_factors,
        forward_factors,
        inverse_factors,
        support_family,
        transport.endpoint_transport_matrix,
        transport.from_certificate,
        transport.to_certificate,
        link_identity,
    )
    verification.overall_ok ||
        throw(ArgumentError("ECP link step segment $(idx) failed exact replay verification"))
    return (;
        index = idx,
        from_path_point,
        to_path_point,
        delta,
        from_column,
        to_column,
        sl2_block,
        sl2_embedding,
        elementary_factors,
        forward_factors,
        inverse_factors,
        support_family,
        endpoint_transport_matrix = transport.endpoint_transport_matrix,
        from_certificate = transport.from_certificate,
        to_certificate = transport.to_certificate,
        link_identity,
        verification,
    )
end

function _ecp_link_step_identity(witness::ECPLinkWitnessRecord, idx::Int, from_column, to_column, delta)
    R = witness.ring
    tail = witness.tail_reductions[idx]
    bezout = witness.bezout_coefficients[idx]
    resultant_value = witness.resultants[idx]
    coverage_multiplier = witness.coverage_multipliers[idx]
    from_path_point = witness.path_points[idx]
    evaluated_tail_coefficients = tuple((
        _ecp_evaluate_polynomial_at_selected_variable(coefficient, witness.selected_variable_index, from_path_point, R)
        for coefficient in tail.lifted_tail_coefficients
    )...)
    evaluated_tilde_G = zero(R)
    for tail_idx in eachindex(evaluated_tail_coefficients)
        evaluated_tilde_G += evaluated_tail_coefficients[tail_idx] * from_column[tail_idx + 1]
    end
    evaluated_f = _ecp_evaluate_polynomial_at_selected_variable(bezout.f, witness.selected_variable_index, from_path_point, R)
    evaluated_h = _ecp_evaluate_polynomial_at_selected_variable(bezout.h, witness.selected_variable_index, from_path_point, R)
    evaluated_resultant = _ecp_evaluate_polynomial_at_selected_variable(resultant_value, witness.selected_variable_index, from_path_point, R)
    bezout_total = evaluated_f * from_column[1] + evaluated_h * evaluated_tilde_G
    divided_differences, divisibility_ok = _ecp_link_step_divided_differences(from_column, to_column, delta, R)
    expected_delta = witness.selected_variable * resultant_value * coverage_multiplier
    return (;
        resultant = resultant_value,
        coverage_multiplier,
        delta,
        expected_delta,
        delta_ok = delta == expected_delta,
        evaluated_tail_coefficients,
        evaluated_tilde_G,
        evaluated_f,
        evaluated_h,
        evaluated_resultant,
        bezout_total,
        bezout_ok = bezout_total == evaluated_resultant,
        divided_differences,
        divisibility_ok,
        overall_ok = delta == expected_delta && bezout_total == evaluated_resultant && divisibility_ok,
    )
end

function _ecp_evaluate_polynomial_at_selected_variable(value, selected_variable_index::Int, point, R)
    substitutions = collect(gens(R))
    substitutions[selected_variable_index] = _coerce_into_ring(R, point, "path point")
    return _coerce_into_ring(R, evaluate(value, substitutions), "evaluated link identity entry")
end

function _ecp_link_step_divided_differences(from_column, to_column, delta, R)
    differences = [to_column[idx] - from_column[idx] for idx in eachindex(from_column)]
    if delta == zero(R)
        return tuple((zero(R) for _ in differences)...), all(iszero, differences)
    end

    quotients = Any[]
    for difference in differences
        quotient = try
            divexact(difference, delta)
        catch err
            err isa InterruptException && rethrow()
            return tuple((zero(R) for _ in differences)...), false
        end
        quotient * delta == difference || return tuple((zero(R) for _ in differences)...), false
        push!(quotients, quotient)
    end
    return tuple(quotients...), true
end

function _ecp_link_step_supported_family(witness::ECPLinkWitnessRecord)
    if _ecp_link_step_matches_gf2_fixture(witness) || _ecp_link_step_matches_qq_fixture(witness)
        return :supplied_fixture_identity_sl2_endpoint_transport
    end
    probe_ids = tuple((probe.id for probe in witness.residue_probes)...)
    throw(ArgumentError("unsupported ECP link step family for supplied link witness probes $(probe_ids)"))
end

function _ecp_link_step_base_characteristic_is(R, value::Integer)
    char = try
        characteristic(base_ring(R))
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    return char == value
end

function _ecp_link_step_probe_ideal_matches(probe, generators)
    return hasproperty(probe, :kind) &&
        hasproperty(probe, :maximal_ideal_generators) &&
        probe.kind == :deterministic_fixture &&
        tuple(probe.maximal_ideal_generators...) == generators
end

function _ecp_link_step_tail_matches(tail, G, lifted_tail_coefficients)
    return hasproperty(tail, :G) &&
        hasproperty(tail, :lifted_tail_coefficients) &&
        hasproperty(tail, :tilde_G) &&
        tail.G == G &&
        tuple(tail.lifted_tail_coefficients...) == lifted_tail_coefficients &&
        tail.tilde_G == G
end

function _ecp_link_step_bezout_matches(bezout, f, h)
    return hasproperty(bezout, :f) &&
        hasproperty(bezout, :h) &&
        bezout.f == f &&
        bezout.h == h
end

function _ecp_link_step_matches_gf2_fixture(witness::ECPLinkWitnessRecord)
    R = witness.ring
    !_ecp_link_step_base_characteristic_is(R, 2) && return false
    ring_generators = gens(R)
    length(ring_generators) == 2 || return false
    x, y = ring_generators[1], ring_generators[2]
    length(witness.original_column) == 3 || return false
    length(witness.residue_probes) == 1 || return false
    length(witness.tail_reductions) == 1 || return false
    length(witness.bezout_coefficients) == 1 || return false

    expected_column = (x + y^2, x * y + x + one(R), x^2 + x * y + y + one(R))
    expected_tail = y * expected_column[2] + expected_column[3]
    return witness.selected_variable_index == 1 &&
        witness.selected_variable == x &&
        witness.selected_monic_index == 1 &&
        witness.original_column == expected_column &&
        _ecp_link_step_probe_ideal_matches(witness.residue_probes[1], (y,)) &&
        _ecp_link_step_tail_matches(witness.tail_reductions[1], expected_tail, (y, one(R))) &&
        witness.resultants == (one(R),) &&
        _ecp_link_step_bezout_matches(witness.bezout_coefficients[1], x, one(R)) &&
        witness.coverage_multipliers == (one(R),) &&
        witness.path_points == (zero(R), x)
end

function _ecp_default_public_link_witness(column, R, selected_variable)
    _ecp_link_step_base_characteristic_is(R, 2) ||
        throw(ArgumentError("unsupported ECP staged public column pipeline route"))
    ring_generators = gens(R)
    length(ring_generators) == 2 ||
        throw(ArgumentError("unsupported ECP staged public column pipeline route"))
    x, y = ring_generators
    selected_variable == x ||
        throw(ArgumentError("unsupported ECP staged public column pipeline route"))
    expected = (x + y^2, x * y + x + one(R), x^2 + x * y + y + one(R))
    tuple(column...) == expected ||
        throw(ArgumentError("unsupported ECP staged public column pipeline route"))
    G = y * column[2] + column[3]
    return (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((; probe_id = :gf2_fixture_probe, G, lifted_tail_coefficients = (y, one(R)), tilde_G = G),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
end

function _ecp_link_step_matches_qq_fixture(witness::ECPLinkWitnessRecord)
    R = witness.ring
    !_ecp_link_step_base_characteristic_is(R, 0) && return false
    ring_generators = gens(R)
    length(ring_generators) == 2 || return false
    x, y = ring_generators[1], ring_generators[2]
    length(witness.original_column) == 3 || return false
    length(witness.residue_probes) == 2 || return false
    length(witness.tail_reductions) == 2 || return false
    length(witness.bezout_coefficients) == 2 || return false

    expected_column = (x^2 + y + one(R), x, y)
    return witness.selected_variable_index == 1 &&
        witness.selected_variable == x &&
        witness.selected_monic_index == 1 &&
        witness.original_column == expected_column &&
        _ecp_link_step_probe_ideal_matches(witness.residue_probes[1], (y,)) &&
        _ecp_link_step_probe_ideal_matches(witness.residue_probes[2], (x,)) &&
        _ecp_link_step_tail_matches(witness.tail_reductions[1], y, (zero(R), one(R))) &&
        _ecp_link_step_tail_matches(witness.tail_reductions[2], x, (one(R), zero(R))) &&
        witness.resultants == (y^2, y + one(R)) &&
        _ecp_link_step_bezout_matches(witness.bezout_coefficients[1], zero(R), y) &&
        _ecp_link_step_bezout_matches(witness.bezout_coefficients[2], one(R), -x) &&
        witness.coverage_multipliers == (one(R), one(R) - y) &&
        witness.path_points == (zero(R), y^2 * x, x)
end

function _ecp_link_step_identity_transport(R, from_column, to_column, link_identity, support_family::Symbol)
    link_identity.overall_ok ||
        throw(ArgumentError("ECP link step requires a replayed link identity"))
    support_family == :supplied_fixture_identity_sl2_endpoint_transport ||
        throw(ArgumentError("unsupported ECP link step family $(support_family)"))
    from_certificate = try
        ecp_column_reduction_certificate(collect(from_column), R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("unsupported ECP link step source path column"))
    end
    to_certificate = try
        ecp_column_reduction_certificate(collect(to_column), R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("unsupported ECP link step target path column"))
    end
    factors = vcat(_ecp_inverse_factor_sequence(to_certificate.factors), from_certificate.factors)
    endpoint_transport_matrix = _factor_sequence_product(factors, R, length(from_column))
    endpoint_transport_matrix * matrix(R, length(from_column), 1, collect(from_column)) ==
        matrix(R, length(to_column), 1, collect(to_column)) ||
        throw(ErrorException("internal error: ECP link step endpoint transport does not map path endpoints"))
    return (;
        from_certificate,
        to_certificate,
        endpoint_transport_matrix,
        factors,
    )
end

function _ecp_default_public_normality_witness(lower_factors, n::Int, R)
    lower_product = _factor_sequence_product(lower_factors, R, n)
    return (;
        source = :supplied_normality_witness,
        conjugator = inv(lower_product),
        sl2_indices = (n, 1),
        sl2_entry = gens(R)[end] + one(R),
    )
end

function _ecp_inverse_elementary_factor(factor)
    R = base_ring(factor)
    n = nrows(factor)
    n == ncols(factor) || throw(ArgumentError("expected a square elementary factor"))
    identity = identity_matrix(R, n)
    positions = [(row, col) for row in 1:n, col in 1:n if factor[row, col] != identity[row, col]]
    length(positions) == 1 || throw(ArgumentError("expected an elementary factor with one off-diagonal entry"))
    row, col = only(positions)
    row != col || throw(ArgumentError("expected an off-diagonal elementary factor"))
    return elementary_matrix(n, row, col, -factor[row, col], R)
end

function _ecp_inverse_factor_sequence(factors)
    return [_ecp_inverse_elementary_factor(factor) for factor in reverse(factors)]
end

function _ecp_link_step_segment_verification(
    R,
    from_column,
    to_column,
    sl2_block,
    sl2_embedding,
    elementary_factors,
    forward_factors,
    inverse_factors,
    support_family,
    endpoint_transport_matrix,
    from_certificate,
    to_certificate,
    link_identity,
)
    n = length(from_column)
    from_matrix = matrix(R, n, 1, collect(from_column))
    to_matrix = matrix(R, n, 1, collect(to_column))
    sl2_block_ok = nrows(sl2_block) == 2 &&
        ncols(sl2_block) == 2 &&
        _same_base_ring(base_ring(sl2_block), R) &&
        det(sl2_block) == one(R) &&
        sl2_block == identity_matrix(R, 2)
    sl2_embedding_ok = sl2_embedding == block_embedding(sl2_block, n, (1, 2))
    support_family_ok = support_family == :supplied_fixture_identity_sl2_endpoint_transport
    endpoint_reductions_ok = verify_ecp_column_reduction(from_certificate) &&
        verify_ecp_column_reduction(to_certificate) &&
        from_certificate.original_column == collect(from_column) &&
        to_certificate.original_column == collect(to_column)
    endpoint_transport_ok = endpoint_transport_matrix == _factor_sequence_product(forward_factors, R, n)
    elementary_factors_ok = _ecp_factor_sequences_equal(elementary_factors, forward_factors)
    forward_map_ok = _apply_reduction_factors(forward_factors, collect(from_column), R) == to_matrix
    inverse_map_ok = _apply_reduction_factors(inverse_factors, collect(to_column), R) == from_matrix
    inverse_sequence_ok = _ecp_factor_sequences_equal(inverse_factors, _ecp_inverse_factor_sequence(forward_factors))
    overall_ok = sl2_block_ok &&
        sl2_embedding_ok &&
        support_family_ok &&
        endpoint_reductions_ok &&
        endpoint_transport_ok &&
        elementary_factors_ok &&
        link_identity.overall_ok &&
        forward_map_ok &&
        inverse_map_ok &&
        inverse_sequence_ok
    return (;
        overall_ok,
        sl2_block_ok,
        sl2_embedding_ok,
        support_family_ok,
        endpoint_reductions_ok,
        endpoint_transport_ok,
        elementary_factors_ok,
        link_identity_ok = link_identity.overall_ok,
        forward_map_ok,
        inverse_map_ok,
        inverse_sequence_ok,
    )
end

function _ecp_link_step_forward_factors(segments)
    factors = Any[]
    for segment in reverse(segments)
        append!(factors, segment.forward_factors)
    end
    return factors
end

function _ecp_link_step_reduction_factors(segments)
    factors = Any[]
    for segment in segments
        append!(factors, segment.inverse_factors)
    end
    return factors
end

function _ecp_link_step_replay_summary(certificate)
    R = certificate.ring
    witness_ok = verify_ecp_link_witness(certificate.link_witness)
    input_ok = certificate.original_column == certificate.link_witness.original_column
    path_points_ok = certificate.path_points == certificate.link_witness.path_points
    recomputed_path_columns = witness_ok ? _ecp_link_step_path_columns(certificate.link_witness) : ()
    path_columns_ok = certificate.path_columns == recomputed_path_columns
    transformed_column_ok = !isempty(recomputed_path_columns) &&
        certificate.transformed_column == certificate.original_column &&
        certificate.transformed_column == recomputed_path_columns[end]
    lower_variable_column_ok = !isempty(recomputed_path_columns) &&
        certificate.lower_variable_column == recomputed_path_columns[1]

    recomputed_segments = try
        path_columns_ok ? _ecp_link_step_segments(certificate.link_witness, recomputed_path_columns) : ()
    catch err
        err isa InterruptException && rethrow()
        ()
    end
    segments_ok = _ecp_link_step_segments_equivalent(certificate.segments, recomputed_segments) &&
        all(segment -> segment.verification.overall_ok, certificate.segments)
    recomputed_forward_factors = _ecp_link_step_forward_factors(certificate.segments)
    recomputed_reduction_factors = _ecp_link_step_reduction_factors(certificate.segments)
    forward_factors_ok = _ecp_factor_sequences_equal(certificate.forward_factors, recomputed_forward_factors)
    reduction_factors_ok = _ecp_factor_sequences_equal(certificate.reduction_factors, recomputed_reduction_factors)

    n = length(certificate.original_column)
    composed_forward_ok = lower_variable_column_ok &&
        transformed_column_ok &&
        _apply_reduction_factors(certificate.forward_factors, collect(certificate.lower_variable_column), R) ==
            matrix(R, n, 1, collect(certificate.transformed_column))
    composed_reduction_ok = lower_variable_column_ok &&
        transformed_column_ok &&
        _apply_reduction_factors(certificate.reduction_factors, collect(certificate.transformed_column), R) ==
            matrix(R, n, 1, collect(certificate.lower_variable_column))
    overall_ok = witness_ok &&
        input_ok &&
        path_points_ok &&
        path_columns_ok &&
        transformed_column_ok &&
        lower_variable_column_ok &&
        segments_ok &&
        forward_factors_ok &&
        reduction_factors_ok &&
        composed_forward_ok &&
        composed_reduction_ok
    return (;
        overall_ok,
        witness_ok,
        input_ok,
        path_points_ok,
        path_columns_ok,
        transformed_column_ok,
        lower_variable_column_ok,
        segments_ok,
        forward_factors_ok,
        reduction_factors_ok,
        composed_forward_ok,
        composed_reduction_ok,
        segment_verifications = tuple((segment.verification for segment in certificate.segments)...),
    )
end

function _ecp_link_step_segments_equivalent(left, right)
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _ecp_link_step_segment_equivalent(left[idx], right[idx]) || return false
    end
    return true
end

function _ecp_link_step_segment_equivalent(left, right)
    return left.index == right.index &&
        left.from_path_point == right.from_path_point &&
        left.to_path_point == right.to_path_point &&
        left.delta == right.delta &&
        left.from_column == right.from_column &&
        left.to_column == right.to_column &&
        left.sl2_block == right.sl2_block &&
        left.sl2_embedding == right.sl2_embedding &&
        _ecp_factor_sequences_equal(left.elementary_factors, right.elementary_factors) &&
        _ecp_factor_sequences_equal(left.forward_factors, right.forward_factors) &&
        _ecp_factor_sequences_equal(left.inverse_factors, right.inverse_factors) &&
        left.support_family == right.support_family &&
        left.endpoint_transport_matrix == right.endpoint_transport_matrix &&
        _ecp_column_certificates_equivalent(left.from_certificate, right.from_certificate) &&
        _ecp_column_certificates_equivalent(left.to_certificate, right.to_certificate) &&
        left.link_identity == right.link_identity &&
        left.verification == right.verification
end

function _ecp_column_certificates_equivalent(left, right)
    return verify_ecp_column_reduction(left) &&
        verify_ecp_column_reduction(right) &&
        left.original_column == right.original_column &&
        left.ring == right.ring &&
        _ecp_factor_sequences_equal(left.factors, right.factors) &&
        left.final_column == right.final_column &&
        left.verification == right.verification
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

function _ecp_staged_column_reduction_replay_summary(certificate)
    route = (:validation, :monicity_forcing, :link_witness, :link_step, :induction_normality)
    original_column_ok = certificate.original_column == certificate.link_step.original_column &&
        certificate.original_column == certificate.induction_normality.original_column
    link_step_ok = verify_ecp_link_step_certificate(certificate.link_step)
    link_witness_ok = link_step_ok && verify_ecp_link_witness(certificate.link_step.link_witness)
    induction_normality_ok = verify_ecp_induction_normality_certificate(certificate.induction_normality)
    monicity_ok = link_witness_ok &&
        certificate.monicity == _ecp_staged_replayed_monicity(certificate.link_step.link_witness)
    lower_reduction_ok = _ecp_staged_lower_reduction_matches(
        certificate.lower_reduction,
        certificate.induction_normality,
    )
    normality_witness_ok =
        certificate.normality_witness == certificate.induction_normality.normality_witness
    factors_match_ok = _ecp_factor_sequences_equal(certificate.factors, certificate.induction_normality.final_factors)
    replayed_final_column = _apply_reduction_factors(
        certificate.factors,
        collect(certificate.original_column),
        certificate.ring,
    )
    final_column_ok = certificate.final_column == replayed_final_column
    target_ok = certificate.final_column == _target_reduced_column(certificate.ring, length(certificate.original_column))
    overall_ok = all((
        original_column_ok,
        link_witness_ok,
        link_step_ok,
        induction_normality_ok,
        monicity_ok,
        lower_reduction_ok,
        normality_witness_ok,
        factors_match_ok,
        final_column_ok,
        target_ok,
    ))
    return (;
        route,
        overall_ok,
        original_column_ok,
        link_witness_ok,
        link_step_ok,
        induction_normality_ok,
        monicity_ok,
        lower_reduction_ok,
        normality_witness_ok,
        factors_match_ok,
        final_column_ok,
        target_ok,
        replayed_final_column,
    )
end

function _ecp_staged_replayed_monicity(link_witness::ECPLinkWitnessRecord)
    return (;
        source = :link_witness,
        variable_order = link_witness.variable_order,
        selected_variable_index = link_witness.selected_variable_index,
        selected_variable = link_witness.selected_variable,
        selected_monic_index = link_witness.selected_monic_index,
        selected_monic_entry = link_witness.selected_monic_entry,
        selected_monic_ok = _ecp_link_witness_replay_summary(link_witness).selected_monic_ok,
    )
end

function _ecp_staged_lower_reduction_matches(lower_reduction, induction_normality)
    induction_normality.lower_reduction_certificate !== nothing &&
        return lower_reduction == induction_normality.lower_reduction_certificate
    _ecp_is_concrete_factor_sequence(lower_reduction) || return false
    return _ecp_factor_sequences_equal(
        collect(lower_reduction),
        induction_normality.lower_variable_factors,
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
        invalid_replay = (; ok = false, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
        required_keys = (
            :kind,
            :input_column,
            :variable_order,
            :variable_index,
            :source_variable_index,
            :source_variable,
            :last_variable_index,
            :target_variable_index,
            :target_variable,
            :shift_power,
            :shift_sign,
            :shift_polynomial,
            :forward_values,
            :inverse_values,
            :forward_substitution,
            :inverse_substitution,
            :transformed_column,
            :selected_monic_index,
            :selected_monic_entry,
            :first_coordinate_strategy,
            :first_coordinate_move_factors,
            :first_coordinate_column,
            :transformed_stage,
            :transformed_factors,
            :inverse_substituted_factors,
            :factors,
            :output_column,
            :variable_change_verification,
        )
        _ecp_stage_keys_ok(stage, required_keys) || return invalid_replay

        ring_gens = collect(gens(R))
        source_variable_index = stage.source_variable_index
        target_variable_index = stage.target_variable_index
        source_index_ok = source_variable_index isa Integer && 1 <= source_variable_index <= length(ring_gens)
        target_index_ok = target_variable_index isa Integer && 1 <= target_variable_index <= length(ring_gens)
        source_index_ok && target_index_ok || return invalid_replay

        variable_order = try
            collect(stage.variable_order)
        catch err
            err isa InterruptException && rethrow()
            return invalid_replay
        end
        isempty(variable_order) && return invalid_replay
        normalized_variable_order = try
            _ecp_normalize_variable_order(R, variable_order)
        catch err
            err isa ArgumentError || rethrow()
            nothing
        end
        normalized_variable_order === nothing && return invalid_replay

        source_variable = stage.source_variable
        target_variable = stage.target_variable
        source_variable == ring_gens[source_variable_index] || return invalid_replay
        target_variable == ring_gens[target_variable_index] || return invalid_replay
        stage.shift_power isa Integer && stage.shift_power >= 0 || return invalid_replay

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
        source_order_index = findfirst(variable -> variable == source_variable, variable_order)
        target_order_index = findfirst(variable -> variable == target_variable, variable_order)
        ok = stage.input_column == _ecp_column_tuple(input_column) &&
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

function _ecp_public_staged_reduction_certificate(
    column,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    supplied_link_witness = nothing,
    lower_reduction = nothing,
    normality_witness = nothing,
)
    _is_laurent_polynomial_ring(R) && return nothing
    _has_at_least_two_generators(R) || return nothing
    normalized_order = _ecp_normalize_variable_order(R, variable_order)
    isempty(normalized_order) && return nothing
    selected_variable = selected_variable === nothing ? first(normalized_order) : selected_variable
    link_witness = if supplied_link_witness === nothing
        try
            _ecp_default_public_link_witness(column, R, selected_variable)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            return nothing
        end
    else
        supplied_link_witness
    end
    return ecp_staged_column_reduction_certificate(
        column,
        R;
        variable_order,
        selected_variable,
        supplied_link_witness = link_witness,
        lower_reduction = lower_reduction,
        normality_witness = normality_witness,
    )
end

function _diagnose_laurent_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return (; supported = true, stage = :unit_entry)

    push!(attempted, :laurent_unit_creation)
    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return (; supported = true, stage = :laurent_unit_creation)

    push!(attempted, :laurent_witness_unit)
    witness_unit = _reduce_via_laurent_witness_unit_certificate(column, R)
    witness_unit !== nothing && return (; supported = true, stage = :laurent_witness_unit)

    push!(attempted, :laurent_normalization)
    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    return _diagnose_polynomial_unimodular_column_reduction(poly_column, P, attempted)
end

function _diagnose_polynomial_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    small = _diagnose_exact_small_column_reduction(column, R, attempted)
    small.supported && return small

    if length(column) > 3
        push!(attempted, :three_entry_block)
        block = _reduce_via_supported_three_block_certificate(column, R)
        block !== nothing && return (; supported = true, stage = :three_entry_block)
    end

    return (; supported = false, stage = nothing)
end

function _diagnose_exact_small_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    supported = _diagnose_supported_unimodular_column_reduction(column, R, attempted)
    supported.supported && return supported

    if _has_at_least_two_generators(R)
        push!(attempted, :monicity_normalization)
        normalized = try
            _reduce_after_monicity_normalization_certificate(column, R)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        normalized !== nothing && return (; supported = true, stage = :monicity_normalization)
    end

    return (; supported = false, stage = nothing)
end

function _diagnose_supported_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return (; supported = true, stage = :unit_entry)

    push!(attempted, :witness_unit)
    witness = try
        _unimodular_witness(column, R)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    witness === nothing && return (; supported = false, stage = nothing)
    witness_unit_idx = findfirst(is_unit, witness)
    witness_unit_idx !== nothing && return (; supported = true, stage = :witness_unit)

    return (; supported = false, stage = nothing)
end

function _column_reduction_error_message(err)
    return sprint(showerror, err)
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

function _unsupported_unimodular_column_reduction_message(column::AbstractVector, R)
    profile = _is_laurent_polynomial_ring(R) ? "Laurent-normalized" : "ordinary polynomial"
    n = length(column)
    return "unsupported exact unimodular column reduction for $(profile) column of length $(n): no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies"
end

function _throw_unsupported_unimodular_column_reduction(column::AbstractVector, R)
    throw(ArgumentError(_unsupported_unimodular_column_reduction_message(column, R)))
end
