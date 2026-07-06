struct ECPColumnReductionCertificate
    original_column
    ring
    stages
    factors::Vector
    final_column
    verification
end

struct ECPInputContext
    column
    ring
    ring_profile
    variables
    variable_order
    column_length::Int
    unimodularity_witness
    selected_variable_index
    selected_variable
    support_classification
    staged_failure_reason
    staged_diagnostic
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

struct ECPMonicityNormalization
    original_column
    ring
    context::ECPInputContext
    variable_order
    selected_variable_index::Int
    selected_variable
    source_variable_index
    source_variable
    shift_power::Int
    shift_sign
    shift_polynomial
    forward_substitution
    inverse_substitution
    transformed_column
    selected_monic_index::Int
    selected_monic_entry
    coordinate_move_factors::Vector
    normalized_column
    inverse_substituted_coordinate_move_factors::Vector
    inverse_substituted_reduction_factors::Vector
    factors::Vector
    verification
end

struct ECPMonicityNormalizationFailure
    context::ECPInputContext
    variable_order
    selected_variable_index::Int
    selected_variable
    selected_monic_index_hint
    max_shift_power::Int
    shift_signs
    search_failure
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

struct ECPLinkWitnessExtractionFailure
    kind::Symbol
    original_column
    ring
    variable_order
    selected_variable_index::Int
    selected_variable
    selected_monic_index::Int
    selected_monic_entry
    max_tail_coefficient_degree::Int
    max_tail_terms::Int
    max_cover_witnesses::Int
    attempted_tail_reductions::Int
    valid_resultants
    message::String
end

struct ECPLinkStepCertificate
    original_column
    ring
    link_witness::ECPLinkWitnessRecord
    path_points
    path_columns
    route_mode::Symbol
    segments
    lower_variable_column
    transformed_column
    forward_factors::Vector
    reduction_factors::Vector
    verification
end

function ECPLinkStepCertificate(
    original_column,
    ring,
    link_witness,
    path_points,
    path_columns,
    segments,
    lower_variable_column,
    transformed_column,
    forward_factors::Vector,
    reduction_factors::Vector,
    verification,
)
    return ECPLinkStepCertificate(
        original_column,
        ring,
        link_witness,
        path_points,
        path_columns,
        :auto,
        segments,
        lower_variable_column,
        transformed_column,
        forward_factors,
        reduction_factors,
        verification,
    )
end

struct ECPInductionNormalityCertificate
    original_column
    ring
    link_step::ECPLinkStepCertificate
    lower_variable_column
    descent_measure
    lower_reduction_certificate
    lower_variable_factors::Vector
    lifted_lower_variable_factors::Vector
    normality_witness
    normality_certificate
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
    return _ecp_column_reduction_certificate_validated(column, R).factors
end

function ecp_column_reduction_certificate(v::AbstractVector, R)
    column = _validated_unimodular_column(v, R)
    return _ecp_column_reduction_certificate_validated(column, R)
end

function diagnose_unimodular_column_reduction(
    v::AbstractVector,
    R;
    allow_general_ecp_pipeline::Bool = true,
    assume_unimodular::Bool = false,
    laurent_large_support_diagnostic_decline::Bool = false,
)
    ring_profile = _column_reduction_ring_profile(R)
    column_length = length(v)
    validation = _diagnose_unimodular_column_preconditions(
        v,
        R,
        ring_profile,
        column_length;
        check_unimodularity = !assume_unimodular,
    )
    validation.status == :ok || return validation.diagnostic
    return _diagnose_unimodular_column_reduction_validated(
        validation.column,
        R,
        ring_profile,
        column_length,
        allow_general_ecp_pipeline = allow_general_ecp_pipeline,
        laurent_large_support_diagnostic_decline = laurent_large_support_diagnostic_decline,
    )
end

function ecp_input_context(
    v::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    unimodularity_witness = nothing,
)
    return ECPInputContext(
        v,
        R;
        variable_order,
        selected_variable,
        unimodularity_witness,
    )
end

function ECPInputContext(
    v::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    unimodularity_witness = nothing,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP input contexts support ordinary polynomial columns only"))

    column = _validated_unimodular_column(v, R)
    ring_profile = _column_reduction_ring_profile(R)
    variables = tuple(gens(R)...)
    normalized_order = tuple(_ecp_normalize_variable_order(R, variable_order)...)
    selected_variable_index, stored_selected_variable =
        _ecp_input_context_selected_variable(R, normalized_order, selected_variable)
    unimodularity_witness === nothing ||
        _check_ecp_input_context_witness_hint(unimodularity_witness, column, R)
    canonical_witness = _unimodular_witness(column, R)
    staged_diagnostic = diagnose_unimodular_column_reduction(
        column,
        R;
        allow_general_ecp_pipeline = false,
    )
    support_classification = staged_diagnostic.status
    staged_failure_reason = _ecp_input_context_staged_failure_reason(staged_diagnostic)

    provisional = ECPInputContext(
        column,
        R,
        ring_profile,
        variables,
        normalized_order,
        length(column),
        canonical_witness,
        selected_variable_index,
        stored_selected_variable,
        support_classification,
        staged_failure_reason,
        staged_diagnostic,
        nothing,
    )
    verification = _ecp_input_context_replay_summary(provisional)
    verification.overall_ok ||
        error("internal ECP input context verification failed")
    context = ECPInputContext(
        provisional.column,
        provisional.ring,
        provisional.ring_profile,
        provisional.variables,
        provisional.variable_order,
        provisional.column_length,
        provisional.unimodularity_witness,
        provisional.selected_variable_index,
        provisional.selected_variable,
        provisional.support_classification,
        provisional.staged_failure_reason,
        provisional.staged_diagnostic,
        verification,
    )
    verify_ecp_input_context(context) ||
        error("internal ECP input context storage verification failed")
    return context
end

function ecp_monicity_normalization(
    context::ECPInputContext;
    selected_variable = nothing,
    selected_monic_index = nothing,
    max_shift_power::Integer = 3,
    shift_signs = (one(context.ring), -one(context.ring)),
)
    verify_ecp_input_context(context) ||
        throw(ArgumentError("ECP monicity normalization requires a verified ECP input context"))
    max_shift_power < 0 && throw(ArgumentError("max_shift_power must be nonnegative"))

    R = context.ring
    resolved_selected_variable = selected_variable === nothing ? context.selected_variable : selected_variable
    resolved_selected_variable === nothing &&
        throw(ArgumentError("ECP monicity normalization requires a selected variable"))
    selected_variable_index = _ecp_selected_variable_index(R, resolved_selected_variable)
    variable_order = _ecp_target_last_variable_order(R, context.variable_order, resolved_selected_variable)
    selected_monic_index_hint = _ecp_selected_monic_index_hint(selected_monic_index)

    direct = _ecp_monicity_normalization_summary(
        collect(context.column),
        R,
        variable_order,
        selected_variable_index;
        selected_monic_index_hint,
        use_full_reduction = true,
    )
    if direct !== nothing
        return _ecp_monicity_normalization_record(context, variable_order, direct)
    end

    signs = tuple((_coerce_into_ring(R, sign, "shift sign") for sign in shift_signs)...)
    attempted = 0
    for source_variable in variable_order[1:(end - 1)]
        source_variable_index = _ecp_generator_index(R, source_variable)
        for shift_power in 1:Int(max_shift_power), shift_sign in signs
            attempted += 1
            summary = _ecp_monicity_normalization_summary(
                collect(context.column),
                R,
                variable_order,
                selected_variable_index;
                source_variable_index,
                shift_power,
                shift_sign,
                selected_monic_index_hint,
                use_full_reduction = true,
            )
            summary === nothing && continue
            return _ecp_monicity_normalization_record(context, variable_order, summary)
        end
    end

    search_failure = _ecp_monicity_search_failure(
        collect(context.column),
        R,
        variable_order,
        max_shift_power,
        signs,
        attempted,
    )
    return ECPMonicityNormalizationFailure(
        context,
        variable_order,
        selected_variable_index,
        gens(R)[selected_variable_index],
        selected_monic_index_hint,
        Int(max_shift_power),
        signs,
        search_failure,
    )
end

function _ecp_input_context_selected_variable(R, variable_order, selected_variable)
    selected_variable === nothing && return nothing, nothing

    selected_variable_index = _ecp_selected_variable_index(R, selected_variable)
    stored_selected_variable = gens(R)[selected_variable_index]
    count(==(stored_selected_variable), variable_order) == 1 ||
        throw(ArgumentError("selected_variable must appear in variable_order"))
    return selected_variable_index, stored_selected_variable
end

function _check_ecp_input_context_witness_hint(witness, column::AbstractVector, R)
    try
        Base.require_one_based_indexing(witness)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("unimodularity_witness must use one-based indexing"))
    end

    length(witness) == length(column) ||
        throw(ArgumentError("unimodularity_witness must have the same length as v"))
    coerced = [
        _coerce_into_ring(R, witness[idx], "unimodularity_witness[$idx]")
        for idx in 1:length(column)
    ]
    _ecp_input_context_witness_total(coerced, column, R) == one(R) ||
        throw(ArgumentError("unimodularity_witness must certify the input column"))
    return coerced
end

function _ecp_input_context_witness_total(witness, column::AbstractVector, R)
    total = zero(R)
    for idx in 1:length(column)
        total += witness[idx] * column[idx]
    end
    return total
end

function _ecp_input_context_staged_failure_reason(diagnostic)
    return diagnostic.status == :unsupported ? diagnostic.failure_code : nothing
end

function _column_reduction_diagnostic(
    status::Symbol,
    failure_code,
    ring_profile,
    column_length::Int,
    attempted_stages,
    message::AbstractString,
    stage_details = (),
)
    return (;
        status,
        failure_code,
        ring_profile,
        column_length,
        attempted_stages = tuple(attempted_stages...),
        message = String(message),
        stage_details = tuple(stage_details...),
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

function _column_reduction_ring_kind(R)
    return _column_reduction_ring_profile(R).kind
end

function _column_reduction_stage_detail(stage::Symbol, R, outcome::Symbol; kwargs...)
    return (;
        stage,
        ring_kind = _column_reduction_ring_kind(R),
        outcome,
        kwargs...,
    )
end

const _LAURENT_DIAGNOSTIC_SUPPORT_TERM_LIMIT = 1000
const _LAURENT_DESCENT_MEASURE_COMPONENTS = (
    :whole_support_count,
    :max_entry_terms,
    :valuation_span,
    :leading_exponent,
    :leading_entry_index,
)
const _LAURENT_DESCENT_OPERATION_FIELDS = (
    :family,
    :target_index,
    :source_index,
    :coefficient,
    :exponent,
    :ring_generators,
)
const _LAURENT_DESCENT_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :operation,
    :before_measure,
    :after_measure,
    :status,
    :replay_status,
    :measure_relation,
)
const _LAURENT_DESCENT_REQUIRED_MEASURE_FIELDS = (
    :status,
    :order,
    :components,
    _LAURENT_DESCENT_MEASURE_COMPONENTS...,
)

function _laurent_descent_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function _require_two_generator_laurent_ring(R; label::AbstractString = "ring")
    _is_laurent_polynomial_ring(R) ||
        throw(ArgumentError("$label must be a Laurent polynomial ring"))
    ngens(R) == 2 ||
        throw(ArgumentError("$label must have exactly two generators"))
    return R
end

function _laurent_descent_ring_generators(R)
    _require_two_generator_laurent_ring(R)
    return Tuple(string.(gens(R)))
end

function _laurent_descent_exponent_tuple(exponent)::Tuple{Int, Int}
    (exponent isa Tuple || exponent isa AbstractVector) ||
        throw(ArgumentError("exponent must be a length-2 tuple or vector"))
    length(exponent) == 2 ||
        throw(ArgumentError("exponent must have exactly two entries"))
    exponent[1] isa Integer ||
        throw(ArgumentError("first exponent must be an integer"))
    exponent[2] isa Integer ||
        throw(ArgumentError("second exponent must be an integer"))
    return (Int(exponent[1]), Int(exponent[2]))
end

function _laurent_descent_entry_support(entry)::Tuple
    iszero(entry) && return ()
    support = Tuple{Int, Int}[]
    for exponent in exponents(entry)
        push!(support, _laurent_descent_exponent_tuple(exponent))
    end
    sort!(support)
    return Tuple(support)
end

function _laurent_descent_support_bounds(support)
    isempty(support) && return nothing
    return (;
        min_exponents = (
            minimum(term[1] for term in support),
            minimum(term[2] for term in support),
        ),
        max_exponents = (
            maximum(term[1] for term in support),
            maximum(term[2] for term in support),
        ),
    )
end

function _laurent_descent_leading_candidate_lt(left, right)::Bool
    left.leading_exponent == right.leading_exponent &&
        return left.entry_index < right.entry_index
    return isless(right.leading_exponent, left.leading_exponent)
end

function _laurent_descent_measure_from_column(column, R; case_id = nothing)
    _require_two_generator_laurent_ring(R)
    generator_names = _laurent_descent_ring_generators(R)
    coerced_column = [
        _coerce_into_ring(R, column[idx], "column[$idx]")
        for idx in 1:length(column)
    ]
    supports = [_laurent_descent_entry_support(entry) for entry in coerced_column]

    whole_support = Set{Tuple{Int, Int}}()
    for support in supports
        union!(whole_support, support)
    end
    whole_bounds = _laurent_descent_support_bounds(whole_support)

    leading_candidates = NamedTuple[]
    for (idx, support) in enumerate(supports)
        isempty(support) && continue
        push!(
            leading_candidates,
            (;
                entry_index = idx,
                leading_exponent = maximum(support),
            ),
        )
    end
    isempty(leading_candidates) &&
        throw(ArgumentError("cannot measure a column with no nonzero entries"))
    ordered_candidates = sort(leading_candidates; lt = _laurent_descent_leading_candidate_lt)
    leading = first(ordered_candidates)
    valuation_span = whole_bounds === nothing ? (0, 0) : (
        whole_bounds.max_exponents[1] - whole_bounds.min_exponents[1],
        whole_bounds.max_exponents[2] - whole_bounds.min_exponents[2],
    )

    measure = (;
        dimension = length(coerced_column),
        ring_generators = generator_names,
        status = :measure_contract,
        order = :lexicographic_minimize,
        components = _LAURENT_DESCENT_MEASURE_COMPONENTS,
        whole_support_count = length(whole_support),
        max_entry_terms = maximum(length, supports; init = 0),
        valuation_span,
        leading_exponent = leading.leading_exponent,
        leading_entry_index = leading.entry_index,
    )
    return case_id === nothing ? measure : merge((; case_id), measure)
end

function _strictly_decreases_laurent_measure(before, after)::Bool
    _laurent_descent_has_fields(before, _LAURENT_DESCENT_REQUIRED_MEASURE_FIELDS) ||
        return false
    _laurent_descent_has_fields(after, _LAURENT_DESCENT_REQUIRED_MEASURE_FIELDS) ||
        return false
    before.status == :measure_contract || return false
    after.status == :measure_contract || return false
    before.order == :lexicographic_minimize || return false
    after.order == before.order || return false
    before.components == _LAURENT_DESCENT_MEASURE_COMPONENTS || return false
    after.components == before.components || return false
    return isless(
        Tuple(getproperty(after, component) for component in _LAURENT_DESCENT_MEASURE_COMPONENTS),
        Tuple(getproperty(before, component) for component in _LAURENT_DESCENT_MEASURE_COMPONENTS),
    )
end

function _laurent_descent_checked_entry_index(index, n::Int, name::AbstractString)::Int
    index isa Integer ||
        throw(ArgumentError("$name must be an integer entry index"))
    checked = Int(index)
    1 <= checked <= n ||
        throw(ArgumentError("$name must be between 1 and $n"))
    return checked
end

function _laurent_descent_operation_status(operation, n::Int, R)::Symbol
    _require_two_generator_laurent_ring(R)
    _laurent_descent_has_fields(operation, _LAURENT_DESCENT_OPERATION_FIELDS) ||
        return :malformed_operation
    operation.family == :entry_addition || return :malformed_operation
    operation.ring_generators == _laurent_descent_ring_generators(R) ||
        return :wrong_ring_generators
    try
        target = _laurent_descent_checked_entry_index(
            operation.target_index,
            n,
            "target_index",
        )
        source = _laurent_descent_checked_entry_index(
            operation.source_index,
            n,
            "source_index",
        )
        target != source || return :malformed_operation
        _laurent_descent_exponent_tuple(operation.exponent)
        _coerce_into_ring(R, operation.coefficient, "coefficient")
    catch err
        err isa InterruptException && rethrow()
        return :malformed_operation
    end
    return :ok
end

function _replay_laurent_elementary_entry_addition(column, R, operation)::Vector
    _require_two_generator_laurent_ring(R)
    operation_status = _laurent_descent_operation_status(operation, length(column), R)
    operation_status == :ok ||
        throw(ArgumentError("invalid Laurent descent operation: $(operation_status)"))

    transformed = [
        _coerce_into_ring(R, column[idx], "column[$idx]")
        for idx in 1:length(column)
    ]
    target = Int(operation.target_index)
    source = Int(operation.source_index)
    exponent = _laurent_descent_exponent_tuple(operation.exponent)
    coefficient = _coerce_into_ring(R, operation.coefficient, "coefficient")
    monomial = coefficient
    for idx in 1:2
        monomial *= gens(R)[idx]^exponent[idx]
    end
    transformed[target] = transformed[target] + monomial * transformed[source]
    return transformed
end

function _validate_laurent_descent_step_certificate(cert, column, R)::Symbol
    try
        _require_two_generator_laurent_ring(R)
        _laurent_descent_has_fields(cert, _LAURENT_DESCENT_CERTIFICATE_FIELDS) ||
            return :missing_certificate_fields
        cert.status == :descent_step_certificate || return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.measure_relation == :strict_decrease || return :wrong_measure_relation
        cert.dimension == length(column) || return :wrong_dimension

        ring_generators = _laurent_descent_ring_generators(R)
        cert.ring_generators == ring_generators || return :wrong_ring_generators

        operation_status = _laurent_descent_operation_status(
            cert.operation,
            length(column),
            R,
        )
        operation_status == :ok || return operation_status

        expected_before_measure = _laurent_descent_measure_from_column(
            column,
            R;
            case_id = cert.case_id,
        )
        cert.before_measure == expected_before_measure ||
            return :stale_before_measure

        after_column = _replay_laurent_elementary_entry_addition(
            column,
            R,
            cert.operation,
        )
        expected_after_measure = _laurent_descent_measure_from_column(
            after_column,
            R;
            case_id = cert.case_id,
        )
        cert.after_measure == expected_after_measure ||
            return :stale_after_measure

        _strictly_decreases_laurent_measure(
            expected_before_measure,
            expected_after_measure,
        ) || return :not_strict_decrease
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :operation_replay_failed
    end
end

function _column_reduction_entry_term_count(entry)
    try
        iszero(entry) && return 0
        return length(collect(coefficients(entry)))
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _column_reduction_max_entry_term_count(column::AbstractVector)
    max_terms = 0
    for entry in column
        term_count = _column_reduction_entry_term_count(entry)
        term_count === nothing && return nothing
        max_terms = max(max_terms, term_count)
    end
    return max_terms
end

function _laurent_diagnostic_large_support_decline(column::AbstractVector)
    max_entry_term_count = _column_reduction_max_entry_term_count(column)
    max_entry_term_count === nothing && return nothing
    max_entry_term_count <= _LAURENT_DIAGNOSTIC_SUPPORT_TERM_LIMIT && return nothing
    return (;
        max_entry_term_count,
        support_term_limit = _LAURENT_DIAGNOSTIC_SUPPORT_TERM_LIMIT,
    )
end

function _laurent_native_ecp_boundary_stage_detail(R)
    return _column_reduction_stage_detail(
        :laurent_native_ecp_boundary,
        R,
        :staged_boundary;
        boundary = :laurent_native_ecp,
        requires_descent_measure = true,
        requires_link_witness = true,
        requires_endpoint_reduction = true,
        requires_laurent_normality_replay = true,
        requires_recursive_peel_integration = true,
        fallback_policy = :diagnostic_only,
    )
end

function _diagnose_unimodular_column_preconditions(
    v::AbstractVector,
    R,
    ring_profile,
    column_length::Int;
    check_unimodularity::Bool = true,
)
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

    check_unimodularity || return (; status = :ok, column)

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

function _diagnose_unimodular_column_reduction_validated(
    column::AbstractVector,
    R,
    ring_profile,
    column_length::Int;
    allow_general_ecp_pipeline::Bool = true,
    laurent_large_support_diagnostic_decline::Bool = false,
)
    attempted = Symbol[]
    details = Any[]
    result = _is_laurent_polynomial_ring(R) ?
        _diagnose_laurent_unimodular_column_reduction(
            column,
            R,
            attempted,
            details;
            laurent_large_support_diagnostic_decline,
        ) :
        _diagnose_polynomial_unimodular_column_reduction(
            column,
            R,
            attempted,
            details;
            allow_general_ecp_pipeline,
        )

    if result.supported
        return _column_reduction_diagnostic(
            :supported,
            nothing,
            ring_profile,
            column_length,
            attempted,
            "exact unimodular column reduction is supported by $(result.stage)",
            details,
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
        details,
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
    normalization::ECPMonicityNormalization;
    supplied_link_witness = nothing,
    max_tail_coefficient_degree::Integer = 1,
    max_tail_terms::Integer = 2,
    max_cover_witnesses::Integer = 3,
)
    verify_ecp_monicity_normalization(normalization) ||
        throw(ArgumentError("ECP link witness extraction requires a verified monicity normalization record"))
    return ecp_link_witness(
        collect(normalization.normalized_column),
        normalization.ring;
        variable_order = normalization.variable_order,
        selected_variable = normalization.selected_variable,
        selected_monic_index = 1,
        supplied_link_witness,
        max_tail_coefficient_degree,
        max_tail_terms,
        max_cover_witnesses,
    )
end

function ecp_link_witness(
    v::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    selected_monic_index::Integer = 1,
    supplied_link_witness = nothing,
    max_tail_coefficient_degree::Integer = 1,
    max_tail_terms::Integer = 2,
    max_cover_witnesses::Integer = 3,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP link witnesses currently support ordinary polynomial columns only"))

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
    supplied_link_witness === nothing &&
        return _ecp_extract_link_witness(
            column,
            R,
            normalized_order,
            selected_variable_index,
            selected_monic_index;
            max_tail_coefficient_degree = Int(max_tail_coefficient_degree),
            max_tail_terms = Int(max_tail_terms),
            max_cover_witnesses = Int(max_cover_witnesses),
        )

    metadata = (; source = _ecp_link_field(supplied_link_witness, :source))
    metadata.source == :supplied_link_witness ||
        throw(ArgumentError("supplied ECP link witness metadata must use source = :supplied_link_witness"))
    try
        return _ecp_link_witness_record_from_data(
            column,
            R,
            normalized_order,
            selected_variable_index,
            selected_monic_index,
            _ecp_link_field(supplied_link_witness, :residue_probes),
            _ecp_link_field(supplied_link_witness, :tail_reductions),
            _ecp_link_field(supplied_link_witness, :resultants),
            _ecp_link_field(supplied_link_witness, :bezout_coefficients),
            _ecp_link_field(supplied_link_witness, :coverage_multipliers),
            _ecp_link_field(supplied_link_witness, :path_points),
            metadata;
            failure_message = "supplied Park-Woodburn ECP link witness data failed exact replay verification",
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError("supplied Park-Woodburn ECP link witness data failed exact replay verification"))
    end
end

function ecp_link_step_certificate(
    v::AbstractVector,
    R;
    link_witness = nothing,
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    selected_monic_index::Integer = 1,
    supplied_link_witness = nothing,
    route_mode::Symbol = :auto,
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
    resolved_route_mode = _ecp_link_step_resolve_route_mode(witness, route_mode)
    segments = _ecp_link_step_segments(witness, path_columns; route_mode = resolved_route_mode)
    forward_factors = _ecp_link_step_forward_factors(segments)
    reduction_factors = _ecp_link_step_reduction_factors(segments)
    provisional = ECPLinkStepCertificate(
        tuple(column...),
        R,
        witness,
        witness.path_points,
        path_columns,
        resolved_route_mode,
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
        provisional.route_mode,
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
    normality_certificate = nothing,
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
    descent_measure = _ecp_induction_descent_measure(link_step, R)
    _ecp_descent_measure_strict(descent_measure) ||
        throw(ArgumentError("ECP induction/normality staged failure: same-context recursive lower-variable call did not strictly reduce selected-variable profile"))
    lower_certificate, lower_factors = try
        _ecp_verified_lower_reduction(lower_reduction, lower_column, R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("ECP induction/normality staged failure: missing lower-variable reduction: $(sprint(showerror, err))"))
    end
    lifted_lower_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in lower_factors]
    resolved_normality_witness = normality_witness === nothing ?
        _ecp_construct_normality_witness(
            lifted_lower_factors,
            length(lower_column),
            R,
            descent_measure.selected_variable,
        ) :
        normality_witness
    normality_rewrite = try
        _ecp_induction_normality_rewrite(
            resolved_normality_witness,
            normality_certificate,
            lower_column,
            lifted_lower_factors,
            R,
        )
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("ECP induction/normality staged failure: missing normality rewrite: $(sprint(showerror, err))"))
    end
    final_factors = _ecp_induction_final_factors(
        lifted_lower_factors,
        normality_rewrite.rewrite_factors,
        link_step.reduction_factors,
        R,
        length(column),
    )
    final_column = _apply_reduction_factors(final_factors, column, R)
    provisional = ECPInductionNormalityCertificate(
        tuple(column...),
        R,
        link_step,
        tuple(lower_column...),
        descent_measure,
        lower_certificate,
        lower_factors,
        lifted_lower_factors,
        resolved_normality_witness,
        normality_rewrite.normality_certificate,
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
        provisional.descent_measure,
        provisional.lower_reduction_certificate,
        provisional.lower_variable_factors,
        provisional.lifted_lower_variable_factors,
        provisional.normality_witness,
        provisional.normality_certificate,
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
    supported = _reduce_supported_unimodular_column_certificate(column, R)
    supported !== nothing && return supported

    if length(column) > 3
        block_factors = _reduce_via_supported_three_block_certificate(column, R)
        block_factors !== nothing && return block_factors

        general_error = nothing
        general = try
            _reduce_via_general_ecp_pipeline_certificate(column, R)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            general_error = err
            nothing
        end
        general !== nothing && return general

        normalization = _reduce_after_monicity_normalization_certificate(column, R)
        normalization !== nothing && return normalization
        general_error === nothing || throw(general_error)
        return nothing # COV_EXCL_LINE
    end

    general_error = nothing
    general = try
        _reduce_via_general_ecp_pipeline_certificate(column, R)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        general_error = err
        nothing
    end
    general !== nothing && return general

    small = _reduce_after_monicity_normalization_certificate(column, R)
    small !== nothing && return small
    general_error === nothing || throw(general_error)

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

function _ecp_general_pipeline_route_metadata(
    link_step::ECPLinkStepCertificate,
    normalized_column_length::Int,
)
    return (;
        source = :general_ecp_pipeline,
        route = :general_ecp_pipeline,
        link_route_mode = link_step.route_mode,
        normalized_column_length,
        segment_support_families = tuple((segment.support_family for segment in link_step.segments)...),
    )
end

function _ecp_substitute_factor_sequence(factors, substitution_map, R)
    substitution_values = collect(_ecp_substitution_map_values(substitution_map))
    return [
        _substitute_matrix_entries(factor, substitution_values, R)
        for factor in factors
    ]
end

function _ecp_link_endpoint_reduction_certificate(column::AbstractVector, R)
    result = _reduce_supported_unimodular_column_certificate(column, R)
    if result === nothing && length(column) > 3
        result = _reduce_via_supported_three_block_certificate(column, R)
    end
    result === nothing && (result = _reduce_after_monicity_normalization_certificate(column, R))
    result !== nothing ||
        throw(ArgumentError("unsupported ECP link step path column endpoint"))
    return _ecp_certificate_from_stage(column, R, result.stage)
end

function _unsupported_general_ecp_pipeline_message(
    column::AbstractVector,
    R,
    detail::AbstractString,
)
    return _unsupported_unimodular_column_reduction_message(
        column,
        R;
        detail = "general ECP pipeline staged failure: $(detail)",
    )
end

function _reduce_via_general_ecp_pipeline_certificate(column::AbstractVector, R)
    _is_laurent_polynomial_ring(R) && return nothing

    context = ecp_input_context(column, R)
    selected_variable = context.selected_variable === nothing ?
        first(context.variable_order) :
        context.selected_variable
    normalization = ecp_monicity_normalization(context; selected_variable)
    if normalization isa ECPMonicityNormalizationFailure
        throw(ArgumentError(_unsupported_general_ecp_pipeline_message(
            column,
            R,
            normalization.search_failure.message,
        )))
    end

    link_witness = ecp_link_witness(normalization)
    if link_witness isa ECPLinkWitnessExtractionFailure
        throw(ArgumentError(_unsupported_general_ecp_pipeline_message( # COV_EXCL_LINE
            column,
            R,
            link_witness.message,
        )))
    end

    normalized_column = collect(normalization.normalized_column)
    link_step = ecp_link_step_certificate(normalized_column, R; link_witness)
    induction = ecp_induction_normality_certificate(normalized_column, R; link_step)

    inverse_substituted_induction_factors = _ecp_substitute_factor_sequence(
        induction.final_factors,
        normalization.inverse_substitution,
        R,
    )
    factors = _checked_reduction_factors(
        vcat(
            inverse_substituted_induction_factors,
            normalization.inverse_substituted_coordinate_move_factors,
        ),
        column,
        R,
        "general ECP pipeline reduction",
    )
    output_column = _apply_reduction_factors(factors, column, R)
    stage = (;
        kind = :ecp_pipeline,
        input_column = _ecp_column_tuple(column),
        route_metadata = _ecp_general_pipeline_route_metadata(
            link_step,
            length(normalization.normalized_column),
        ),
        context,
        normalization,
        link_witness,
        link_step,
        induction_normality = induction,
        inverse_substituted_induction_factors,
        factors,
        output_column,
    )
    return (; factors, stage)
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

function _ecp_link_witness_record_from_data(
    column,
    R,
    normalized_order,
    selected_variable_index::Int,
    selected_monic_index::Int,
    residue_probes,
    tail_reductions,
    resultants,
    bezout_coefficients,
    coverage_multipliers,
    path_points,
    metadata;
    failure_message = "Park-Woodburn ECP link witness data failed exact replay verification",
)
    replay_record = ECPLinkWitnessRecord(
        tuple(column...),
        R,
        tuple(normalized_order...),
        selected_variable_index,
        gens(R)[selected_variable_index],
        selected_monic_index,
        column[selected_monic_index],
        tuple(residue_probes...),
        tuple(tail_reductions...),
        tuple(resultants...),
        tuple(bezout_coefficients...),
        tuple(coverage_multipliers...),
        tuple(path_points...),
        metadata,
        nothing,
    )
    verification = _ecp_link_witness_replay_summary(replay_record)
    verification.overall_ok || throw(ArgumentError(failure_message))
    stored = ECPLinkWitnessRecord(
        replay_record.original_column,
        replay_record.ring,
        replay_record.variable_order,
        replay_record.selected_variable_index,
        replay_record.selected_variable,
        replay_record.selected_monic_index,
        replay_record.selected_monic_entry,
        replay_record.residue_probes,
        replay_record.tail_reductions,
        replay_record.resultants,
        replay_record.bezout_coefficients,
        replay_record.coverage_multipliers,
        replay_record.path_points,
        replay_record.metadata,
        verification,
    )
    verify_ecp_link_witness(stored) ||
        throw(ArgumentError("stored Park-Woodburn ECP link witness data failed exact replay verification"))
    return stored
end

function _ecp_link_small_scalars(R)
    scalars = Any[one(R)]
    negative_one = -one(R)
    negative_one == scalars[1] || push!(scalars, negative_one)
    return tuple(scalars...)
end

function _ecp_link_exponent_tuples(var_count::Int, max_degree::Int)
    var_count >= 0 || throw(ArgumentError("var_count must be nonnegative"))
    max_degree >= 0 || throw(ArgumentError("max_degree must be nonnegative"))

    exponents = NTuple{var_count, Int}[]
    current = Vector{Int}(undef, var_count)
    function build(idx::Int, remaining::Int)
        if idx > var_count
            push!(exponents, tuple(current...))
            return
        end
        for exponent in 0:remaining
            current[idx] = exponent
            build(idx + 1, remaining - exponent)
        end
    end
    build(1, max_degree)
    sort!(exponents; by = entry -> (sum(entry), entry))
    return tuple(exponents...)
end

function _ecp_link_monomial_basis(R, max_degree::Int)
    max_degree >= 0 || throw(ArgumentError("max_degree must be nonnegative"))

    variables = gens(R)
    basis = Any[]
    for exponents in _ecp_link_exponent_tuples(length(variables), max_degree)
        monomial = one(R)
        for idx in eachindex(variables)
            exponent = exponents[idx]
            exponent == 0 && continue
            monomial *= variables[idx]^exponent
        end
        push!(basis, monomial)
    end
    return tuple(basis...)
end

function _ecp_link_combinations(indices, width::Int)
    width >= 0 || throw(ArgumentError("width must be nonnegative"))
    collected = collect(indices)
    width == 0 && return ((),)
    width > length(collected) && return ()

    combinations = Tuple[]
    current = Vector{Any}(undef, width)
    function build(start_idx::Int, depth::Int)
        if depth > width
            push!(combinations, tuple(current...))
            return
        end
        last_start = length(collected) - (width - depth)
        for idx in start_idx:last_start
            current[depth] = collected[idx]
            build(idx + 1, depth + 1)
        end
    end
    build(1, 1)
    return tuple(combinations...)
end

function _ecp_link_tail_reduction_candidates(
    tail_entries,
    R;
    max_tail_coefficient_degree::Int,
    max_tail_terms::Int,
)
    max_tail_coefficient_degree >= 0 ||
        throw(ArgumentError("max_tail_coefficient_degree must be nonnegative"))
    max_tail_terms >= 0 || throw(ArgumentError("max_tail_terms must be nonnegative"))

    tail_count = length(tail_entries)
    tail_count == 0 && return ()

    atoms = Any[]
    for monomial in _ecp_link_monomial_basis(R, max_tail_coefficient_degree)
        for scalar in _ecp_link_small_scalars(R)
            atom = scalar * monomial
            atom == zero(R) && continue
            any(existing -> existing == atom, atoms) || push!(atoms, atom)
        end
    end

    candidates = NamedTuple[]
    seen_coefficients = Tuple[]
    max_width = min(max_tail_terms, tail_count)
    for width in 1:max_width
        for support in _ecp_link_combinations(1:tail_count, width)
            chosen = Vector{Any}(undef, width)
            function build_atom_assignment(depth::Int)
                if depth > width
                    lifted_tail_coefficients = ntuple(tail_idx -> begin
                        support_position = findfirst(==(tail_idx), support)
                        support_position === nothing ? zero(R) : chosen[support_position]
                    end, tail_count)
                    any(existing -> existing == lifted_tail_coefficients, seen_coefficients) && return
                    G = zero(R)
                    for tail_idx in 1:tail_count
                        G += lifted_tail_coefficients[tail_idx] * tail_entries[tail_idx]
                    end
                    G == zero(R) && return
                    push!(seen_coefficients, lifted_tail_coefficients)
                    push!(candidates, (; lifted_tail_coefficients, G))
                    return
                end
                for atom in atoms
                    chosen[depth] = atom
                    build_atom_assignment(depth + 1)
                end
            end
            build_atom_assignment(1)
        end
    end
    return tuple(candidates...)
end

function _ecp_link_coordinates_tuple(coordinates_value, R, expected_length::Int, label::AbstractString)
    expected_length >= 0 || throw(ArgumentError("expected_length must be nonnegative"))
    nrows(coordinates_value) == 1 ||
        throw(ArgumentError("$(label) coordinates must have exactly one row"))
    ncols(coordinates_value) == expected_length ||
        throw(ArgumentError("$(label) coordinates must have length $(expected_length)"))
    return tuple((
        _coerce_into_ring(R, coordinates_value[1, idx], "$(label) coordinate[$idx]")
        for idx in 1:expected_length
    )...)
end

function _ecp_link_bezout_for_resultant(v1, G, resultant_value, R)
    bezout_ideal = ideal(R, [v1, G])
    resultant_value in bezout_ideal ||
        throw(ArgumentError("resultant must belong to the ideal generated by v1 and G"))
    coordinates_value = Oscar.coordinates(resultant_value, bezout_ideal)
    f, h = _ecp_link_coordinates_tuple(coordinates_value, R, 2, "link Bezout")
    f * v1 + h * G == resultant_value ||
        throw(ArgumentError("link Bezout coordinates failed exact verification"))
    return (; f, h)
end

function _ecp_link_cover_multipliers(resultants, R)
    cover_ideal = ideal(R, collect(resultants))
    one(R) in cover_ideal ||
        throw(ArgumentError("resultants do not generate the unit ideal"))
    coordinates_value = Oscar.coordinates(one(R), cover_ideal)
    multipliers = _ecp_link_coordinates_tuple(coordinates_value, R, length(resultants), "link cover")
    coverage_total = zero(R)
    for idx in eachindex(multipliers)
        coverage_total += multipliers[idx] * resultants[idx]
    end
    coverage_total == one(R) ||
        throw(ArgumentError("link cover coordinates failed exact verification"))
    return multipliers
end

function _ecp_extract_link_witness(
    column,
    R,
    normalized_order,
    selected_variable_index::Int,
    selected_monic_index::Int;
    max_tail_coefficient_degree::Int,
    max_tail_terms::Int,
    max_cover_witnesses::Int,
)
    max_tail_coefficient_degree >= 0 ||
        throw(ArgumentError("max_tail_coefficient_degree must be nonnegative"))
    max_tail_terms >= 0 || throw(ArgumentError("max_tail_terms must be nonnegative"))
    max_cover_witnesses >= 0 || throw(ArgumentError("max_cover_witnesses must be nonnegative"))
    selected_monic_index == 1 ||
        throw(ArgumentError("Park-Woodburn ECP link witnesses require the selected monic entry to be first"))
    _is_monic_in_variable(column[selected_monic_index], R, selected_variable_index) ||
        throw(ArgumentError("selected first entry must be monic in the selected variable"))

    v1 = column[selected_monic_index]
    tail_entries = column[2:end]
    candidates = _ecp_link_tail_reduction_candidates(
        tail_entries,
        R;
        max_tail_coefficient_degree,
        max_tail_terms,
    )
    attempted_tail_reductions = length(candidates)
    valid_candidates = NamedTuple[]
    for candidate in candidates
        resultant_value = resultant(v1, candidate.G, selected_variable_index)
        resultant_value == zero(R) && continue
        bezout = _ecp_link_bezout_for_resultant(v1, candidate.G, resultant_value, R)
        push!(valid_candidates, merge(candidate, (; resultant = resultant_value, bezout)))
    end

    selected_variable = gens(R)[selected_variable_index]
    max_width = min(max_cover_witnesses, length(valid_candidates))
    for width in 1:max_width
        for cover_indices in _ecp_link_combinations(1:length(valid_candidates), width)
            cover_resultants = tuple((valid_candidates[idx].resultant for idx in cover_indices)...)
            cover_multipliers = try
                _ecp_link_cover_multipliers(cover_resultants, R)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                continue
            end

            residue_probes = Any[]
            tail_reductions = Any[]
            resultants = Any[]
            bezout_coefficients = Any[]
            coverage_multipliers = Any[]
            path_points = Any[zero(R)]
            cumulative_path_point = zero(R)
            maximal_ideal_generators = tuple(normalized_order...)
            for (local_idx, candidate_idx) in enumerate(cover_indices)
                candidate = valid_candidates[candidate_idx]
                probe_id = Symbol("bounded_tail_probe_$(candidate_idx)")
                push!(residue_probes, (;
                    id = probe_id,
                    kind = :bounded_tail_combination,
                    maximal_ideal_generators,
                ))
                push!(tail_reductions, (;
                    probe_id,
                    G = candidate.G,
                    lifted_tail_coefficients = candidate.lifted_tail_coefficients,
                    tilde_G = candidate.G,
                ))
                push!(resultants, candidate.resultant)
                push!(bezout_coefficients, candidate.bezout)
                push!(coverage_multipliers, cover_multipliers[local_idx])
                cumulative_path_point +=
                    cover_multipliers[local_idx] * candidate.resultant * selected_variable
                push!(path_points, cumulative_path_point)
            end

            return _ecp_link_witness_record_from_data(
                column,
                R,
                normalized_order,
                selected_variable_index,
                selected_monic_index,
                residue_probes,
                tail_reductions,
                resultants,
                bezout_coefficients,
                coverage_multipliers,
                path_points,
                (; source = :extracted_link_witness);
                failure_message = "extracted Park-Woodburn ECP link witness data failed exact replay verification",
            )
        end
    end

    return ECPLinkWitnessExtractionFailure(
        :link_witness_cover_not_proved,
        tuple(column...),
        R,
        tuple(normalized_order...),
        selected_variable_index,
        selected_variable,
        selected_monic_index,
        column[selected_monic_index],
        max_tail_coefficient_degree,
        max_tail_terms,
        max_cover_witnesses,
        attempted_tail_reductions,
        tuple((candidate.resultant for candidate in valid_candidates)...),
        "bounded Park-Woodburn ECP link witness extraction did not prove a cover",
    )
end

function _ecp_selected_variable_index(R, variable)
    ring_gens = collect(gens(R))
    idx = _ecp_variable_order_match_index(ring_gens, variable)
    idx === nothing && throw(ArgumentError("selected_variable must be a generator of R"))
    return idx
end

function _ecp_target_last_variable_order(R, variable_order, selected_variable)
    normalized_order = tuple(_ecp_normalize_variable_order(R, variable_order)...)
    selected_variable_index = _ecp_selected_variable_index(R, selected_variable)
    selected_generator = gens(R)[selected_variable_index]
    count(==(selected_generator), normalized_order) == 1 ||
        throw(ArgumentError("selected_variable must appear in variable_order"))
    prefix = [variable for variable in normalized_order if variable != selected_generator]
    push!(prefix, selected_generator)
    return tuple(prefix...)
end

function _ecp_selected_monic_index_hint(selected_monic_index)
    selected_monic_index === nothing && return nothing
    selected_monic_index isa Integer ||
        throw(ArgumentError("selected_monic_index must be an integer or nothing"))
    Int(selected_monic_index) >= 1 ||
        throw(ArgumentError("selected_monic_index must be positive"))
    return Int(selected_monic_index)
end

function _ecp_monicity_normalization_record(context::ECPInputContext, variable_order, summary)
    provisional = ECPMonicityNormalization(
        tuple(context.column...),
        context.ring,
        context,
        tuple(variable_order...),
        summary.target_variable_index,
        summary.target_variable,
        summary.source_variable_index,
        summary.source_variable,
        summary.shift_power,
        summary.shift_sign,
        summary.shift_polynomial,
        summary.forward_substitution,
        summary.inverse_substitution,
        summary.transformed_column,
        summary.selected_monic_index,
        summary.selected_monic_entry,
        summary.coordinate_move_factors,
        summary.normalized_column,
        summary.inverse_substituted_coordinate_move_factors,
        summary.inverse_substituted_reduction_factors,
        summary.inverse_substituted_factors,
        nothing,
    )
    verification = _ecp_monicity_normalization_replay_summary(provisional)
    verification.overall_ok ||
        error("internal ECP monicity normalization verification failed")
    record = ECPMonicityNormalization(
        provisional.original_column,
        provisional.ring,
        provisional.context,
        provisional.variable_order,
        provisional.selected_variable_index,
        provisional.selected_variable,
        provisional.source_variable_index,
        provisional.source_variable,
        provisional.shift_power,
        provisional.shift_sign,
        provisional.shift_polynomial,
        provisional.forward_substitution,
        provisional.inverse_substitution,
        provisional.transformed_column,
        provisional.selected_monic_index,
        provisional.selected_monic_entry,
        provisional.coordinate_move_factors,
        provisional.normalized_column,
        provisional.inverse_substituted_coordinate_move_factors,
        provisional.inverse_substituted_reduction_factors,
        provisional.factors,
        verification,
    )
    verify_ecp_monicity_normalization(record) ||
        error("internal ECP monicity normalization storage verification failed")
    return record
end

function _ecp_resolve_selected_monic_index(
    transformed_column,
    R,
    target_variable_index::Int,
    selected_monic_index_hint,
)
    selected_monic_index_hint === nothing &&
        return _ecp_first_monic_entry_index(transformed_column, R, target_variable_index)
    1 <= selected_monic_index_hint <= length(transformed_column) || return nothing
    _is_monic_in_variable(
        transformed_column[selected_monic_index_hint],
        R,
        target_variable_index,
    ) || return nothing
    return selected_monic_index_hint
end

function _ecp_monicity_normalization_full_reduction_certificate(column::AbstractVector, R)
    length(column) >= 3 || return nothing

    supported = _reduce_supported_unimodular_column_certificate(column, R)
    supported !== nothing && return supported

    if length(column) > 3
        block = _reduce_via_supported_three_block_certificate(column, R)
        block !== nothing && return block
    end

    return _reduce_after_monicity_normalization_certificate(column, R)
end

function _ecp_monicity_normalization_summary(
    column::AbstractVector,
    R,
    variable_order,
    target_variable_index::Int;
    source_variable_index = nothing,
    shift_power::Integer = 0,
    shift_sign = zero(R),
    selected_monic_index_hint = nothing,
    use_full_reduction::Bool = false,
)
    ring_gens = collect(gens(R))
    normalized_order = tuple(_ecp_normalize_variable_order(R, variable_order)...)
    target_variable = ring_gens[target_variable_index]
    normalized_order[end] == target_variable || return nothing

    source_variable = nothing
    forward_values = copy(ring_gens)
    inverse_values = copy(ring_gens)
    shift_power = Int(shift_power)
    shift_power >= 0 || return nothing
    shift_polynomial = zero(R)
    if source_variable_index !== nothing
        source_variable_index isa Integer || return nothing
        1 <= source_variable_index <= length(ring_gens) || return nothing
        source_variable_index == target_variable_index && return nothing
        shift_power >= 1 || return nothing
        source_variable = ring_gens[source_variable_index]
        shift_polynomial = shift_sign * target_variable^shift_power
        forward_values[source_variable_index] = source_variable + shift_polynomial
        inverse_values[source_variable_index] = source_variable - shift_polynomial
    else
        shift_power == 0 || return nothing
        iszero(shift_sign) || return nothing
    end

    transformed = [
        _coerce_into_ring(R, evaluate(entry, forward_values), "substituted column entry")
        for entry in column
    ]
    selected_monic_index = _ecp_resolve_selected_monic_index(
        transformed,
        R,
        target_variable_index,
        selected_monic_index_hint,
    )
    selected_monic_index === nothing && return nothing
    selected_monic_entry = transformed[selected_monic_index]

    coordinate_move_factors = _ecp_first_coordinate_move_factors(
        length(column),
        selected_monic_index,
        R,
    )
    normalized_column = _ecp_matrix_column_to_tuple(
        _apply_reduction_factors(coordinate_move_factors, transformed, R),
    )
    normalized_column[1] == selected_monic_entry || return nothing

    transformed_result = use_full_reduction ?
        _ecp_monicity_normalization_full_reduction_certificate(
            collect(normalized_column),
            R,
        ) :
        _reduce_supported_unimodular_column_certificate(
            collect(normalized_column),
            R,
        )
    transformed_result === nothing && return nothing

    transformed_factors = _checked_reduction_factors(
        vcat(transformed_result.factors, coordinate_move_factors),
        transformed,
        R,
        "monicity normalization transformed reduction",
    )
    transformed_output = _apply_reduction_factors(
        transformed_factors,
        transformed,
        R,
    )
    inverse_substituted_coordinate_move_factors = [
        _substitute_matrix_entries(factor, inverse_values, R)
        for factor in coordinate_move_factors
    ]
    inverse_substituted_reduction_factors = [
        _substitute_matrix_entries(factor, inverse_values, R)
        for factor in transformed_result.factors
    ]
    inverse_substituted_factors = vcat(
        inverse_substituted_reduction_factors,
        inverse_substituted_coordinate_move_factors,
    )
    inverse_substituted_factors = _checked_reduction_factors(
        inverse_substituted_factors,
        column,
        R,
        "monicity normalization reduction",
    )
    original_output = _apply_reduction_factors(
        inverse_substituted_factors,
        column,
        R,
    )
    target = _target_reduced_column(R, length(column))
    forward_substitution = _ecp_substitution_map_tuple(ring_gens, forward_values)
    inverse_substitution = _ecp_substitution_map_tuple(ring_gens, inverse_values)
    variable_change_verification = (;
        selected_monic_ok = _is_monic_in_variable(selected_monic_entry, R, target_variable_index),
        substitution_inverse_ok = _ecp_substitution_maps_are_inverse(
            R,
            forward_substitution,
            inverse_substitution,
        ),
        transformed_reduction_ok = transformed_output == target,
        original_reduction_ok = original_output == target,
    )
    return (;
        variable_order = normalized_order,
        source_variable_index,
        source_variable,
        target_variable_index,
        target_variable,
        shift_power,
        shift_sign,
        shift_polynomial,
        forward_values = tuple(forward_values...),
        inverse_values = tuple(inverse_values...),
        forward_substitution,
        inverse_substitution,
        transformed_column = tuple(transformed...),
        selected_monic_index,
        selected_monic_entry,
        first_coordinate_strategy = _ecp_first_coordinate_strategy(selected_monic_index),
        coordinate_move_factors,
        normalized_column,
        transformed_stage = transformed_result.stage,
        transformed_reduction_factors = transformed_result.factors,
        transformed_factors,
        inverse_substituted_coordinate_move_factors,
        inverse_substituted_reduction_factors,
        inverse_substituted_factors,
        transformed_output,
        original_output,
        variable_change_verification,
    )
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
    summary = _ecp_monicity_normalization_summary(
        column,
        R,
        variable_order,
        target_variable_index;
        source_variable_index,
        shift_power,
        shift_sign,
    )
    summary === nothing && return nothing
    stage = (;
        kind = :monicity_normalization,
        input_column = _ecp_column_tuple(column),
        variable_order = summary.variable_order,
        variable_index = source_variable_index,
        source_variable_index,
        source_variable = summary.source_variable,
        last_variable_index = target_variable_index,
        target_variable_index,
        target_variable = summary.target_variable,
        shift_power,
        shift_sign,
        shift_polynomial = summary.shift_polynomial,
        forward_values = summary.forward_values,
        inverse_values = summary.inverse_values,
        forward_substitution = summary.forward_substitution,
        inverse_substitution = summary.inverse_substitution,
        transformed_column = summary.transformed_column,
        selected_monic_index = summary.selected_monic_index,
        selected_monic_entry = summary.selected_monic_entry,
        first_coordinate_strategy = summary.first_coordinate_strategy,
        first_coordinate_move_factors = summary.coordinate_move_factors,
        first_coordinate_column = summary.normalized_column,
        transformed_stage = summary.transformed_stage,
        transformed_factors = summary.transformed_factors,
        inverse_substituted_coordinate_move_factors = summary.inverse_substituted_coordinate_move_factors,
        inverse_substituted_reduction_factors = summary.inverse_substituted_reduction_factors,
        inverse_substituted_factors = summary.inverse_substituted_factors,
        factors = summary.inverse_substituted_factors,
        output_column = summary.original_output,
        variable_change_verification = summary.variable_change_verification,
    )
    return (; factors = summary.inverse_substituted_factors, stage)
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

function _laurent_row_preconditioning_specs(column::AbstractVector, R)
    _is_laurent_polynomial_ring(R) || return ()
    length(gens(R)) == 2 || return ()

    if length(column) == 16
        return ((;
            target_index = 1,
            source_indices = (10,),
            coefficient_strategy = :fixed_coefficients,
            coefficients = (one(R),),
            max_nonzero_coefficients = 1,
        ),)
    end

    if length(column) == 15
        return ((;
            target_index = 1,
            source_indices = Tuple(2:15),
            coefficient_strategy = :target_unit_laurent_linear_synthesis,
            coefficients = (),
            max_nonzero_coefficients = 14,
        ),)
    end

    return ()
end

function _prefer_laurent_row_preconditioning_before_base(column::AbstractVector, R)::Bool
    _is_laurent_polynomial_ring(R) || return false
    findfirst(is_unit, column) === nothing || return false
    specs = _laurent_row_preconditioning_specs(column, R)
    return any(spec -> spec.coefficient_strategy == :target_unit_laurent_linear_synthesis, specs)
end

function _reduce_laurent_unimodular_column_base_certificate(column::AbstractVector, R)
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

    poly_result = try
        _reduce_polynomial_unimodular_column_exact_certificate(poly_column, P)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        nothing
    end
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

function _laurent_row_preconditioning_solve_failure(err)::Bool
    return _laurent_witness_solve_failure(err)
end

function _laurent_row_preconditioning_fixed_coefficients(spec, R)
    return tuple((
        _coerce_into_ring(R, coeff, "row preconditioning coefficient")
        for coeff in spec.coefficients
    )...)
end

function _laurent_row_preconditioning_synthesis_coefficients(
    column::AbstractVector,
    R,
    target_idx::Int,
    source_indices,
    ;
    solver = solve_laurent_linear,
)
    isempty(source_indices) && return nothing

    A = matrix(R, 1, length(source_indices), [column[idx] for idx in source_indices])
    B = matrix(R, 1, 1, [one(R) - column[target_idx]])
    solution = try
        solver(A, B)
    catch err
        err isa InterruptException && rethrow()
        _laurent_row_preconditioning_solve_failure(err) && return nothing
        rethrow()
    end

    A * solution == B || return nothing
    return tuple((solution[idx, 1] for idx in 1:nrows(solution))...)
end

function _laurent_row_preconditioning_coefficients(column::AbstractVector, R, spec)
    target_idx = Int(spec.target_index)
    source_indices = tuple(spec.source_indices...)
    strategy = spec.coefficient_strategy
    if strategy == :fixed_coefficients
        length(spec.coefficients) == length(source_indices) || return nothing
        return _laurent_row_preconditioning_fixed_coefficients(spec, R)
    elseif strategy == :target_unit_laurent_linear_synthesis
        return _laurent_row_preconditioning_synthesis_coefficients(
            column,
            R,
            target_idx,
            source_indices,
        )
    end
    return nothing
end

function _laurent_row_preconditioning_source_order_ok(source_indices, allowed_source_indices)::Bool
    isempty(source_indices) && return false

    next_allowed = 1
    for source_idx in source_indices
        found = false
        for allowed_idx in next_allowed:length(allowed_source_indices)
            if source_idx == allowed_source_indices[allowed_idx]
                next_allowed = allowed_idx + 1
                found = true
                break
            end
        end
        found || return false
    end
    return true
end

function _laurent_row_preconditioning_target_unit_equation_ok(
    column::AbstractVector,
    R,
    target_idx::Int,
    source_indices,
    coefficients,
)::Bool
    total = zero(R)
    for (source_idx, coeff) in zip(source_indices, coefficients)
        total += coeff * column[source_idx]
    end
    return total == one(R) - column[target_idx]
end

function _laurent_row_preconditioning_stage_spec_ok(
    column::AbstractVector,
    R,
    target_idx::Int,
    source_indices,
    coefficients,
    coefficient_strategy::Symbol,
)::Bool
    length(source_indices) == length(coefficients) || return false
    isempty(source_indices) && return false
    any(coeff -> coeff == zero(R), coefficients) && return false

    for spec in _laurent_row_preconditioning_specs(column, R)
        spec.target_index == target_idx || continue
        spec.coefficient_strategy == coefficient_strategy || continue
        length(source_indices) <= spec.max_nonzero_coefficients || continue

        allowed_source_indices = tuple(spec.source_indices...)
        _laurent_row_preconditioning_source_order_ok(source_indices, allowed_source_indices) ||
            continue

        if coefficient_strategy == :fixed_coefficients
            fixed = _laurent_row_preconditioning_fixed_coefficients(spec, R)
            nonzero_pairs = [
                (source_idx, coeff)
                for (source_idx, coeff) in zip(allowed_source_indices, fixed)
                if coeff != zero(R)
            ]
            tuple((pair[1] for pair in nonzero_pairs)...) == source_indices || continue
            tuple((pair[2] for pair in nonzero_pairs)...) == coefficients || continue
            return true
        elseif coefficient_strategy == :target_unit_laurent_linear_synthesis
            _laurent_row_preconditioning_target_unit_equation_ok(
                column,
                R,
                target_idx,
                source_indices,
                coefficients,
            ) || continue
            return true
        end
    end

    return false
end

function _laurent_row_preconditioning_candidate(column::AbstractVector, R)
    n = length(column)
    for spec in _laurent_row_preconditioning_specs(column, R)
        spec.target_index isa Integer || continue
        target_idx = Int(spec.target_index)
        1 <= target_idx <= n || continue

        source_indices = tuple((Int(source_idx) for source_idx in spec.source_indices)...)
        all(source_idx -> 1 <= source_idx <= n && source_idx != target_idx, source_indices) ||
            continue
        length(unique(source_indices)) == length(source_indices) || continue

        coefficients = _laurent_row_preconditioning_coefficients(column, R, spec)
        coefficients === nothing && continue
        length(coefficients) == length(source_indices) || continue

        nonzero_pairs = [
            (source_idx, coeff)
            for (source_idx, coeff) in zip(source_indices, coefficients)
            if coeff != zero(R)
        ]
        isempty(nonzero_pairs) && continue
        length(nonzero_pairs) <= spec.max_nonzero_coefficients || continue

        accepted_source_indices = tuple((pair[1] for pair in nonzero_pairs)...)
        accepted_coefficients = tuple((pair[2] for pair in nonzero_pairs)...)
        precondition_factors = [
            elementary_matrix(n, target_idx, source_idx, coeff, R)
            for (source_idx, coeff) in nonzero_pairs
        ]
        precondition_factor = _factor_sequence_product(precondition_factors, R, n)
        transformed_column = collect(_ecp_matrix_column_to_tuple(
            precondition_factor * matrix(R, n, 1, collect(column)),
        ))
        transformed_result =
            _reduce_laurent_unimodular_column_base_certificate(transformed_column, R)
        transformed_result === nothing && continue
        transformed_certificate =
            _ecp_certificate_from_stage(transformed_column, R, transformed_result.stage)
        return (;
            target_index = target_idx,
            source_index = accepted_source_indices[1],
            source_indices = accepted_source_indices,
            coefficient = accepted_coefficients[1],
            coefficients = accepted_coefficients,
            coefficient_strategy = spec.coefficient_strategy,
            precondition_factors,
            precondition_factor,
            transformed_column,
            transformed_certificate,
        )
    end

    return nothing
end

function _reduce_via_laurent_elementary_row_preconditioning_certificate(column::AbstractVector, R)
    candidate = _laurent_row_preconditioning_candidate(column, R)
    candidate === nothing && return nothing

    factors = _checked_reduction_factors(
        vcat(candidate.transformed_certificate.factors, candidate.precondition_factors),
        column,
        R,
        "Laurent elementary row-preconditioning reduction",
    )
    stage = (;
        kind = :laurent_elementary_row_preconditioning,
        input_column = _ecp_column_tuple(column),
        target_index = candidate.target_index,
        source_index = candidate.source_index,
        source_indices = candidate.source_indices,
        coefficient = candidate.coefficient,
        coefficients = candidate.coefficients,
        coefficient_strategy = candidate.coefficient_strategy,
        precondition_factors = candidate.precondition_factors,
        precondition_factor = candidate.precondition_factor,
        transformed_column = tuple(candidate.transformed_column...),
        transformed_certificate = candidate.transformed_certificate,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end

function _reduce_laurent_unimodular_column_certificate(column::AbstractVector, R)
    prefer_row_preconditioning = _prefer_laurent_row_preconditioning_before_base(column, R)
    if prefer_row_preconditioning
        row_preconditioned = _reduce_via_laurent_elementary_row_preconditioning_certificate(column, R)
        row_preconditioned !== nothing && return row_preconditioned
    end

    base = _reduce_laurent_unimodular_column_base_certificate(column, R)
    base !== nothing && return base

    if !prefer_row_preconditioning
        row_preconditioned = _reduce_via_laurent_elementary_row_preconditioning_certificate(column, R)
        row_preconditioned !== nothing && return row_preconditioned
    end

    return nothing
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

function verify_ecp_input_context(context)::Bool
    try
        replay = _ecp_input_context_replay_summary(context)
        return replay.overall_ok && context.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function verify_ecp_monicity_normalization(record)::Bool
    try
        replay = _ecp_monicity_normalization_replay_summary(record)
        return replay.overall_ok && record.verification == replay
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
    certificate = ecp_column_reduction_certificate(lower_column, R)
    _ecp_factor_sequences_equal(factors, certificate.factors) ||
        throw(ArgumentError("lower-variable factor sequence must match the verified lower ECP certificate"))
    return certificate, certificate.factors
end

function _ecp_is_concrete_factor_sequence(value)
    value isa AbstractVector && return true
    value isa Tuple && !(value isa NamedTuple) && return true
    return false
end

function _ecp_selected_variable_profile(column, R, selected_variable_index::Int)
    return tuple((Int(degree(entry, selected_variable_index)) for entry in column)...)
end

function _ecp_induction_descent_measure(link_step::ECPLinkStepCertificate, R)
    verify_ecp_link_step_certificate(link_step) ||
        throw(ArgumentError("ECP induction/normality descent requires a verified link-step certificate"))
    _same_base_ring(link_step.ring, R) ||
        throw(ArgumentError("ECP induction/normality descent ring must match the link-step ring"))
    selected_variable_index = link_step.link_witness.selected_variable_index
    selected_variable = gens(R)[selected_variable_index]
    parent_profile = _ecp_selected_variable_profile(link_step.original_column, R, selected_variable_index)
    lower_profile = _ecp_selected_variable_profile(link_step.lower_variable_column, R, selected_variable_index)
    componentwise_nonincreasing = all(
        lower_profile[idx] <= parent_profile[idx]
        for idx in eachindex(parent_profile)
    )
    strict_descent = componentwise_nonincreasing &&
        any(lower_profile[idx] < parent_profile[idx] for idx in eachindex(parent_profile))
    return (;
        variable_count = ngens(R),
        column_length = length(link_step.original_column),
        selected_variable_index,
        selected_variable,
        parent_profile,
        lower_profile,
        componentwise_nonincreasing,
        strict_descent,
    )
end

function _ecp_descent_measure_strict(measure)::Bool
    try
        return measure.componentwise_nonincreasing === true &&
            measure.strict_descent === true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _ecp_construct_normality_witness(lifted_lower_factors, n::Int, R, selected_variable)
    lower_product = _factor_sequence_product(lifted_lower_factors, R, n)
    entry = _coerce_into_ring(R, selected_variable + one(R), "constructed normality witness sl2_entry")
    entry != zero(R) ||
        throw(ArgumentError("constructed normality witness produced an identity SL_2 contribution"))
    return (;
        source = :constructed_normality_witness,
        conjugator = inv(lower_product),
        sl2_indices = (n, 1),
        sl2_entry = entry,
    )
end

function _ecp_induction_final_factors(lifted_lower_factors, rewrite_factors, link_reduction_factors, R, n::Int)
    return vcat(lifted_lower_factors, rewrite_factors, link_reduction_factors)
end

function _ecp_induction_normality_rewrite(normality_witness, lower_column, lifted_lower_factors, R)
    return _ecp_induction_normality_rewrite(normality_witness, nothing, lower_column, lifted_lower_factors, R)
end

function _ecp_induction_normality_rewrite(normality_witness, normality_certificate, lower_column, lifted_lower_factors, R)
    normality_witness === nothing &&
        throw(ArgumentError("ECP induction/normality requires explicit normality witness data"))
    _ecp_normality_witness_keys_ok(normality_witness) ||
        throw(ArgumentError("normality witness must contain source, conjugator, sl2_indices, and sl2_entry"))
    normality_witness.source in (:supplied_normality_witness, :constructed_normality_witness) ||
        throw(ArgumentError("normality witness must use a supported ECP normality witness source"))

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
    nested_certificate = try
        _ecp_induction_normality_certificate_from_inputs(
            normality_certificate,
            conjugator,
            fixed_index,
            moving_index,
            entry,
            R,
            n,
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError("normality witness could not be rewritten into elementary factors"))
    end
    rewrite_factors = nested_certificate.factors
    rewrite_product = _factor_sequence_product(rewrite_factors, R, n)
    expected_rewrite_product = conjugator * sl2_embedding * lower_product
    lower_matrix = matrix(R, n, 1, collect(lower_column))
    fixed_lower_column_ok = rewrite_product * lower_matrix == lower_matrix
    rewrite_product_ok = rewrite_product == expected_rewrite_product
    return (;
        source = normality_witness.source,
        conjugator,
        lower_product,
        sl2_indices,
        fixed_index,
        moving_index,
        sl2_entry = entry,
        sl2_block,
        sl2_embedding,
        normality_certificate = nested_certificate,
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

function _ecp_induction_normality_certificate_from_inputs(
    supplied_certificate,
    conjugator,
    fixed_index::Int,
    moving_index::Int,
    entry,
    R,
    n::Int,
)
    expected = realize_conjugate_elementary_certificate(conjugator, fixed_index, moving_index, entry)
    supplied_certificate === nothing && return expected
    verify_conjugate_elementary_certificate(supplied_certificate) ||
        throw(ArgumentError("normality certificate does not verify"))
    _same_base_ring(supplied_certificate.ring, R) ||
        throw(ArgumentError("normality certificate ring must match the ECP ring"))
    supplied_certificate.n == n ||
        throw(ArgumentError("normality certificate dimension must match the ECP column"))
    _ecp_conjugated_normality_certificates_match(supplied_certificate, expected) ||
        throw(ArgumentError("normality certificate does not match the supplied witness data"))
    return supplied_certificate
end

function _ecp_conjugated_normality_certificates_match(actual, expected)
    return actual.n == expected.n &&
        actual.A == expected.A &&
        actual.i == expected.i &&
        actual.j == expected.j &&
        actual.a == expected.a &&
        _same_base_ring(actual.ring, expected.ring) &&
        actual.determinant == expected.determinant &&
        actual.inverse_A == expected.inverse_A &&
        actual.elementary_matrix == expected.elementary_matrix &&
        actual.conjugation_convention == expected.conjugation_convention &&
        actual.conjugation_target == expected.conjugation_target &&
        actual.v == expected.v &&
        actual.w == expected.w &&
        actual.g == expected.g &&
        actual.rank_one_certificate.n == expected.rank_one_certificate.n &&
        actual.rank_one_certificate.v == expected.rank_one_certificate.v &&
        actual.rank_one_certificate.w == expected.rank_one_certificate.w &&
        actual.rank_one_certificate.g == expected.rank_one_certificate.g &&
        _same_base_ring(actual.rank_one_certificate.ring, expected.rank_one_certificate.ring) &&
        actual.rank_one_certificate.orthogonality == expected.rank_one_certificate.orthogonality &&
        actual.rank_one_certificate.bezout == expected.rank_one_certificate.bezout &&
        actual.rank_one_certificate.cohn_coefficients == expected.rank_one_certificate.cohn_coefficients &&
        actual.rank_one_certificate.factors == expected.rank_one_certificate.factors &&
        actual.rank_one_certificate.target == expected.rank_one_certificate.target &&
        actual.rank_one_certificate.product == expected.rank_one_certificate.product &&
        actual.rank_one_certificate.verification == expected.rank_one_certificate.verification &&
        actual.factors == expected.factors &&
        actual.product == expected.product &&
        actual.verification == expected.verification
end

function _ecp_input_context_replay_summary(context)
    R = context.ring
    one_based_indexing_ok = try
        Base.require_one_based_indexing(context.column)
        true
    catch err
        err isa InterruptException && rethrow()
        false
    end
    replayed_column_length = try
        length(context.column)
    catch err
        err isa InterruptException && rethrow()
        -1
    end
    ordinary_polynomial_ring_ok = try
        !_is_laurent_polynomial_ring(R)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    replayed_ring_profile = _column_reduction_ring_profile(R)
    ring_profile_ok = context.ring_profile == replayed_ring_profile
    replayed_variables = try
        tuple(gens(R)...)
    catch err
        err isa InterruptException && rethrow()
        ()
    end
    variables_ok = context.variables == replayed_variables
    replayed_variable_order = try
        tuple(_ecp_normalize_variable_order(R, context.variable_order)...)
    catch err
        err isa InterruptException && rethrow()
        ()
    end
    variable_order_ok = context.variable_order == replayed_variable_order
    column_length_ok = context.column_length == replayed_column_length && replayed_column_length >= 3

    replayed_column = try
        one_based_indexing_ok && replayed_column_length >= 0 || throw(ArgumentError("invalid context column"))
        [
            _coerce_into_ring(R, context.column[idx], "context.column[$idx]")
            for idx in 1:replayed_column_length
        ]
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    column_ok = replayed_column !== nothing && context.column == replayed_column
    unimodular_ok = try
        replayed_column !== nothing && is_unimodular_column(replayed_column, R)
    catch err
        err isa InterruptException && rethrow()
        false
    end

    replayed_witness = try
        unimodular_ok ? _unimodular_witness(replayed_column, R) : Any[]
    catch err
        err isa InterruptException && rethrow()
        Any[]
    end
    witness_one_based_indexing_ok = try
        Base.require_one_based_indexing(context.unimodularity_witness)
        true
    catch err
        err isa InterruptException && rethrow()
        false
    end
    witness_length_ok = try
        witness_one_based_indexing_ok &&
            length(context.unimodularity_witness) == replayed_column_length
    catch err
        err isa InterruptException && rethrow()
        false
    end
    replayed_stored_witness = try
        witness_length_ok || throw(ArgumentError("invalid context witness"))
        [
            _coerce_into_ring(
                R,
                context.unimodularity_witness[idx],
                "context.unimodularity_witness[$idx]",
            )
            for idx in 1:replayed_column_length
        ]
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    witness_coercion_ok = replayed_stored_witness !== nothing
    witness_identity_ok = try
        witness_coercion_ok &&
            replayed_column !== nothing &&
            _ecp_input_context_witness_total(replayed_stored_witness, replayed_column, R) == one(R)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    canonical_witness_ok = context.unimodularity_witness == replayed_witness
    unimodularity_witness_ok = witness_one_based_indexing_ok &&
        witness_length_ok &&
        witness_coercion_ok &&
        witness_identity_ok &&
        canonical_witness_ok

    selected_variable_index_ok, replayed_selected_variable_index = try
        if context.selected_variable_index === nothing
            true, nothing
        elseif context.selected_variable_index isa Integer
            idx = Int(context.selected_variable_index)
            1 <= idx <= length(replayed_variables), idx
        else
            false, nothing
        end
    catch err
        err isa InterruptException && rethrow()
        false, nothing
    end
    replayed_selected_variable = selected_variable_index_ok &&
        replayed_selected_variable_index !== nothing ?
        replayed_variables[replayed_selected_variable_index] :
        nothing
    selected_variable_generator_index = try
        context.selected_variable === nothing ? nothing :
            _ecp_selected_variable_index(R, context.selected_variable)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    selected_variable_order_ok = try
        context.selected_variable === nothing ||
            (variable_order_ok && count(==(context.selected_variable), context.variable_order) == 1)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    selected_variable_ok = if context.selected_variable_index === nothing ||
            context.selected_variable === nothing
        context.selected_variable_index === nothing && context.selected_variable === nothing
    else
        selected_variable_index_ok &&
            context.selected_variable_index == replayed_selected_variable_index &&
            selected_variable_generator_index == replayed_selected_variable_index &&
            context.selected_variable == replayed_selected_variable &&
            selected_variable_order_ok
    end

    replayed_staged_diagnostic = try
        ordinary_polynomial_ring_ok && column_ok && unimodular_ok ?
            diagnose_unimodular_column_reduction(
                replayed_column,
                R;
                allow_general_ecp_pipeline = false,
            ) :
            nothing
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    replayed_staged_failure_reason = replayed_staged_diagnostic === nothing ?
        nothing :
        _ecp_input_context_staged_failure_reason(replayed_staged_diagnostic)
    support_classification_ok = replayed_staged_diagnostic !== nothing &&
        context.support_classification == replayed_staged_diagnostic.status
    staged_failure_reason_ok = replayed_staged_diagnostic !== nothing &&
        context.staged_failure_reason == replayed_staged_failure_reason
    staged_diagnostic_ok = replayed_staged_diagnostic !== nothing &&
        context.staged_diagnostic == replayed_staged_diagnostic &&
        support_classification_ok &&
        staged_failure_reason_ok

    overall_ok = ordinary_polynomial_ring_ok &&
        one_based_indexing_ok &&
        column_ok &&
        ring_profile_ok &&
        variables_ok &&
        variable_order_ok &&
        column_length_ok &&
        unimodular_ok &&
        unimodularity_witness_ok &&
        selected_variable_ok &&
        staged_diagnostic_ok
    return (;
        overall_ok,
        ordinary_polynomial_ring_ok,
        one_based_indexing_ok,
        column_ok,
        ring_profile_ok,
        variables_ok,
        variable_order_ok,
        column_length_ok,
        unimodular_ok,
        witness_one_based_indexing_ok,
        witness_length_ok,
        witness_coercion_ok,
        witness_identity_ok,
        canonical_witness_ok,
        unimodularity_witness_ok,
        selected_variable_index_ok,
        selected_variable_order_ok,
        selected_variable_ok,
        support_classification_ok,
        staged_failure_reason_ok,
        staged_diagnostic_ok,
        replayed_column,
        replayed_ring_profile,
        replayed_variables,
        replayed_variable_order,
        replayed_column_length,
        replayed_witness,
        replayed_selected_variable_index,
        replayed_selected_variable,
        replayed_staged_failure_reason,
        replayed_staged_diagnostic,
    )
end

function _ecp_monicity_normalization_replay_summary(record)
    context_ok = verify_ecp_input_context(record.context)
    context_column = context_ok ? tuple(record.context.column...) : ()
    context_ring = context_ok ? record.context.ring : nothing
    input_ok = context_ok &&
        record.original_column == context_column &&
        record.ring == context_ring

    selected_variable_index_ok = try
        1 <= record.selected_variable_index <= length(gens(record.ring))
    catch err
        err isa InterruptException && rethrow()
        false
    end
    replayed_selected_variable = selected_variable_index_ok ?
        gens(record.ring)[record.selected_variable_index] :
        nothing
    selected_variable_ok = selected_variable_index_ok &&
        record.selected_variable == replayed_selected_variable

    replayed_variable_order = try
        context_ok && selected_variable_ok ?
            _ecp_target_last_variable_order(
                record.ring,
                record.context.variable_order,
                record.selected_variable,
            ) :
            ()
    catch err
        err isa InterruptException && rethrow()
        ()
    end
    variable_order_ok = record.variable_order == replayed_variable_order

    summary = try
        input_ok && selected_variable_ok && variable_order_ok ?
            _ecp_monicity_normalization_summary(
                collect(record.original_column),
                record.ring,
                record.variable_order,
                record.selected_variable_index;
                source_variable_index = record.source_variable_index,
                shift_power = record.shift_power,
                shift_sign = record.shift_sign,
                selected_monic_index_hint = record.selected_monic_index,
                use_full_reduction = true,
            ) :
            nothing
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    summary_ok = summary !== nothing
    replayed_target = _target_reduced_column(record.ring, length(record.original_column))
    exact_reduction_ok = summary_ok && summary.original_output == replayed_target
    substitution_inverse_ok = summary_ok &&
        summary.variable_change_verification.substitution_inverse_ok

    fields_ok = summary_ok &&
        record.source_variable_index == summary.source_variable_index &&
        record.source_variable == summary.source_variable &&
        record.shift_power == summary.shift_power &&
        record.shift_sign == summary.shift_sign &&
        record.shift_polynomial == summary.shift_polynomial &&
        record.forward_substitution == summary.forward_substitution &&
        record.inverse_substitution == summary.inverse_substitution &&
        record.transformed_column == summary.transformed_column &&
        record.selected_monic_index == summary.selected_monic_index &&
        record.selected_monic_entry == summary.selected_monic_entry &&
        _ecp_factor_sequences_equal(record.coordinate_move_factors, summary.coordinate_move_factors) &&
        record.normalized_column == summary.normalized_column &&
        _ecp_factor_sequences_equal(
            record.inverse_substituted_coordinate_move_factors,
            summary.inverse_substituted_coordinate_move_factors,
        ) &&
        _ecp_factor_sequences_equal(
            record.inverse_substituted_reduction_factors,
            summary.inverse_substituted_reduction_factors,
        ) &&
        _ecp_factor_sequences_equal(record.factors, summary.inverse_substituted_factors)

    overall_ok = context_ok &&
        input_ok &&
        selected_variable_index_ok &&
        selected_variable_ok &&
        variable_order_ok &&
        summary_ok &&
        substitution_inverse_ok &&
        fields_ok &&
        exact_reduction_ok
    return (;
        overall_ok,
        context_ok,
        input_ok,
        selected_variable_index_ok,
        selected_variable_ok,
        variable_order_ok,
        summary_ok,
        substitution_inverse_ok,
        fields_ok,
        exact_reduction_ok,
        replayed_selected_variable,
        replayed_variable_order,
        replayed_target,
    )
end

function _ecp_induction_normality_replay_summary(certificate)
    R = certificate.ring
    n = length(certificate.original_column)
    link_step_ok = verify_ecp_link_step_certificate(certificate.link_step)
    input_ok = link_step_ok && certificate.original_column == certificate.link_step.original_column
    lower_variable_column_ok = link_step_ok &&
        certificate.lower_variable_column == certificate.link_step.lower_variable_column
    descent_measure = try
        link_step_ok ? _ecp_induction_descent_measure(certificate.link_step, R) : nothing
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    descent_measure_ok = descent_measure !== nothing &&
        certificate.descent_measure == descent_measure
    descent_strict_ok = descent_measure !== nothing &&
        _ecp_descent_measure_strict(descent_measure)

    lower_certificate, lower_factors = try
        _ecp_verified_lower_reduction(
            certificate.lower_reduction_certificate,
            collect(certificate.lower_variable_column),
            R,
        )
    catch err
        err isa InterruptException && rethrow()
        nothing, Any[]
    end
    lower_reduction_certificate_ok = lower_certificate !== nothing &&
        certificate.lower_reduction_certificate == lower_certificate &&
        verify_ecp_column_reduction(certificate.lower_reduction_certificate)
    lower_reduction_ok = lower_reduction_certificate_ok &&
        _ecp_factor_sequences_equal(certificate.lower_variable_factors, lower_factors)
    lifted_lower_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in lower_factors]
    lifted_lower_factors_ok = lower_reduction_certificate_ok &&
        _ecp_factor_sequences_equal(certificate.lifted_lower_variable_factors, lifted_lower_factors)

    normality_rewrite = try
        _ecp_induction_normality_rewrite(
            certificate.normality_witness,
            certificate.normality_certificate,
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
    normality_certificate_ok = normality_rewrite !== nothing &&
        certificate.normality_certificate == normality_rewrite.normality_certificate &&
        verify_conjugate_elementary_certificate(certificate.normality_certificate)

    expected_final_factors = normality_rewrite === nothing ?
        Any[] :
        _ecp_induction_final_factors(
            certificate.lifted_lower_variable_factors,
            normality_rewrite.rewrite_factors,
            certificate.link_step.reduction_factors,
            R,
            n,
        )
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
        descent_measure_ok &&
        descent_strict_ok &&
        lower_reduction_certificate_ok &&
        lower_reduction_ok &&
        lifted_lower_factors_ok &&
        normality_rewrite_ok &&
        normality_certificate_ok &&
        final_factors_ok &&
        final_column_ok &&
        final_reduction_ok &&
        final_factors_elementary_ok
    return (;
        overall_ok,
        link_step_ok,
        input_ok,
        lower_variable_column_ok,
        descent_measure_ok,
        descent_strict_ok,
        lower_reduction_certificate_ok,
        lower_reduction_ok,
        lifted_lower_factors_ok,
        normality_rewrite_ok,
        normality_certificate_ok,
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
        # Link-step certificates can store the zero elementary factor E_ij(0).
        isempty(positions) && return true
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

    metadata_ok = hasproperty(record.metadata, :source) &&
        record.metadata.source in (:supplied_link_witness, :extracted_link_witness)
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

function _ecp_link_step_segments(witness::ECPLinkWitnessRecord, path_columns; route_mode::Symbol = :auto)
    resolved_route_mode = _ecp_link_step_resolve_route_mode(witness, route_mode)
    if resolved_route_mode == :direct_elementary
        return (_ecp_link_step_direct_segment(witness, path_columns),)
    end
    return tuple((
        _ecp_link_step_segment(witness, idx, path_columns; route_mode = resolved_route_mode)
        for idx in eachindex(witness.tail_reductions)
    )...)
end

function _ecp_link_step_direct_segment(witness::ECPLinkWitnessRecord, path_columns)
    R = witness.ring
    n = length(witness.original_column)
    from_path_point = first(witness.path_points)
    to_path_point = last(witness.path_points)
    delta = to_path_point - from_path_point
    from_column = first(path_columns)
    to_column = last(path_columns)
    segment_identities = tuple((
        _ecp_link_step_identity(
            witness,
            idx,
            path_columns[idx],
            path_columns[idx + 1],
            witness.path_points[idx + 1] - witness.path_points[idx],
        )
        for idx in eachindex(witness.tail_reductions)
    )...)
    link_identity = (;
        segment_identities,
        overall_ok = all(identity -> identity.overall_ok, segment_identities),
    )
    support_family = :direct_elementary_endpoint_transport
    transport = _ecp_link_step_endpoint_transport(R, from_column, to_column, link_identity, support_family)

    sl2_block = _ecp_link_step_embedded_sl2_block(transport.sl3_route_matrices, R, n, (1, 2))
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
        transport.sl3_route_matrices,
        transport.sl3_route_certificates,
        transport.sl3_route_factor_groups,
        transport.sl3_route_metadata,
    )
    verification.overall_ok ||
        throw(ArgumentError("ECP link step direct segment failed exact replay verification"))
    return (;
        index = 1,
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
        sl3_route_matrices = transport.sl3_route_matrices,
        sl3_route_certificates = transport.sl3_route_certificates,
        sl3_route_factor_groups = transport.sl3_route_factor_groups,
        sl3_route_metadata = transport.sl3_route_metadata,
        link_identity,
        verification,
    )
end

function _ecp_link_step_segment(witness::ECPLinkWitnessRecord, idx::Int, path_columns; route_mode::Symbol = :auto)
    R = witness.ring
    n = length(witness.original_column)
    from_path_point = witness.path_points[idx]
    to_path_point = witness.path_points[idx + 1]
    delta = to_path_point - from_path_point
    from_column = path_columns[idx]
    to_column = path_columns[idx + 1]
    link_identity = _ecp_link_step_identity(witness, idx, from_column, to_column, delta)
    support_family = _ecp_link_step_supported_family(witness, route_mode)
    transport = _ecp_link_step_endpoint_transport(R, from_column, to_column, link_identity, support_family)

    sl2_block = _ecp_link_step_embedded_sl2_block(transport.sl3_route_matrices, R, n, (1, 2))
    sl2_embedding = block_embedding(sl2_block, n, (1, 2))
    sl3_route_metadata = _ecp_link_step_embedded_route_metadata(
        transport.sl3_route_metadata,
        transport.sl3_route_matrices,
        sl2_embedding,
        sl2_block,
        (1, 2),
    )
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
        transport.sl3_route_matrices,
        transport.sl3_route_certificates,
        transport.sl3_route_factor_groups,
        sl3_route_metadata,
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
        sl3_route_matrices = transport.sl3_route_matrices,
        sl3_route_certificates = transport.sl3_route_certificates,
        sl3_route_factor_groups = transport.sl3_route_factor_groups,
        sl3_route_metadata,
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

function _ecp_link_step_resolve_route_mode(witness::ECPLinkWitnessRecord, route_mode::Symbol)
    route_mode in (:auto, :legacy_fixture, :polynomial_sl3, :direct_elementary) ||
        throw(ArgumentError("unsupported ECP link step route_mode $(route_mode)"))
    if route_mode == :auto
        return (
            _ecp_link_step_matches_gf2_fixture(witness) ||
            _ecp_link_step_matches_qq_fixture(witness)
        ) ?
            :legacy_fixture :
            (length(witness.original_column) > 3 ? :direct_elementary : :polynomial_sl3)
    end
    return route_mode
end

function _ecp_link_step_supported_family(witness::ECPLinkWitnessRecord, route_mode::Symbol)
    resolved = _ecp_link_step_resolve_route_mode(witness, route_mode)
    if resolved == :legacy_fixture
        if _ecp_link_step_matches_gf2_fixture(witness) || _ecp_link_step_matches_qq_fixture(witness)
            return :supplied_fixture_identity_sl2_endpoint_transport
        end
        probe_ids = tuple((probe.id for probe in witness.residue_probes)...)
        throw(ArgumentError("unsupported ECP legacy fixture link step family for supplied link witness probes $(probe_ids)"))
    end
    resolved == :direct_elementary && return :direct_elementary_endpoint_transport
    return :polynomial_sl3_route_endpoint_transport
end

function _ecp_link_step_supported_family(witness::ECPLinkWitnessRecord)
    return _ecp_link_step_supported_family(witness, :auto)
end

function _ecp_link_step_extract_embedded_sl2_block(route_matrix, R, n::Int, indices)
    nrows(route_matrix) == n && ncols(route_matrix) == n || return nothing
    _same_base_ring(base_ring(route_matrix), R) || return nothing
    i, j = indices
    block = matrix(R, 2, 2, [route_matrix[i, i], route_matrix[i, j], route_matrix[j, i], route_matrix[j, j]])
    return block_embedding(block, n, indices) == route_matrix ? block : nothing
end

function _ecp_link_step_embedded_sl2_block(route_matrices, R, n::Int, indices)
    identity_block = identity_matrix(R, 2)
    for route_matrix in route_matrices
        block = _ecp_link_step_extract_embedded_sl2_block(route_matrix, R, n, indices)
        block === nothing && continue
        block == identity_block && continue
        return block
    end
    return identity_block
end

function _ecp_link_step_embedded_route_metadata(route_metadata, route_matrices, sl2_embedding, sl2_block, indices)
    length(route_metadata) == length(route_matrices) ||
        throw(ArgumentError("route metadata must align with route matrices"))
    return tuple((
        route_matrices[idx] == sl2_embedding ?
            merge(route_metadata[idx], (;
                embedded_block_indices = indices,
                embedded_block_matrix = sl2_block,
            )) :
            route_metadata[idx]
        for idx in eachindex(route_metadata)
    )...)
end

function _ecp_link_step_route_certificate(factor)
    certificate = _polynomial_factorization_route_certificate(
        factor;
        allow_recursive_column_peel = false,
    )
    certificate.status == :supported ||
        throw(ArgumentError("ECP link step SL_3 route obligation is staged"))
    _verify_polynomial_factorization_route_certificate(certificate) ||
        throw(ArgumentError("ECP link step SL_3 route certificate does not verify"))
    certificate.product == factor ||
        throw(ArgumentError("ECP link step SL_3 route certificate product does not match its obligation"))
    return certificate
end

function _ecp_link_step_route_metadata(route_certificates)
    return tuple((
        (;
            source = :polynomial_factorization_route_certificate,
            obligation_index = idx,
            route = route_certificates[idx].route,
            status = route_certificates[idx].status,
            factor_count = length(route_certificates[idx].factors),
        )
        for idx in eachindex(route_certificates)
    )...)
end

function _ecp_link_step_endpoint_transport(R, from_column, to_column, link_identity, support_family::Symbol)
    link_identity.overall_ok ||
        throw(ArgumentError("ECP link step requires a replayed link identity"))
    from_certificate = try
        _ecp_link_endpoint_reduction_certificate(collect(from_column), R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("unsupported ECP link step source path column"))
    end
    to_certificate = try
        _ecp_link_endpoint_reduction_certificate(collect(to_column), R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("unsupported ECP link step target path column"))
    end
    raw_factors = vcat(_ecp_inverse_factor_sequence(to_certificate.factors), from_certificate.factors)

    if support_family == :supplied_fixture_identity_sl2_endpoint_transport
        endpoint_transport_matrix = _factor_sequence_product(raw_factors, R, length(from_column))
        endpoint_transport_matrix * matrix(R, length(from_column), 1, collect(from_column)) ==
            matrix(R, length(to_column), 1, collect(to_column)) ||
            throw(ErrorException("internal error: ECP link step endpoint transport does not map path endpoints"))
        return (; from_certificate, to_certificate, endpoint_transport_matrix, factors = raw_factors,
            sl3_route_matrices = (), sl3_route_certificates = (), sl3_route_factor_groups = (),
            sl3_route_metadata = ())
    elseif support_family == :direct_elementary_endpoint_transport
        endpoint_transport_matrix = _factor_sequence_product(raw_factors, R, length(from_column))
        endpoint_transport_matrix * matrix(R, length(from_column), 1, collect(from_column)) ==
            matrix(R, length(to_column), 1, collect(to_column)) ||
            throw(ErrorException("internal error: ECP link step endpoint transport does not map path endpoints"))
        return (; from_certificate, to_certificate, endpoint_transport_matrix, factors = raw_factors,
            sl3_route_matrices = (), sl3_route_certificates = (), sl3_route_factor_groups = (),
            sl3_route_metadata = ())
    elseif support_family == :polynomial_sl3_route_endpoint_transport
        route_certificates = Any[]
        for factor in raw_factors
            push!(route_certificates, _ecp_link_step_route_certificate(factor))
        end
        sl3_route_certificates = tuple(route_certificates...)
        sl3_route_matrices = tuple((cert.matrix for cert in sl3_route_certificates)...)
        sl3_route_factor_groups = tuple((tuple(cert.factors...) for cert in sl3_route_certificates)...)
        factors = Any[]
        for group in sl3_route_factor_groups
            append!(factors, group)
        end
        endpoint_transport_matrix = _factor_sequence_product(factors, R, length(from_column))
        endpoint_transport_matrix * matrix(R, length(from_column), 1, collect(from_column)) ==
            matrix(R, length(to_column), 1, collect(to_column)) ||
            throw(ErrorException("internal error: ECP link step endpoint transport does not map path endpoints"))
        sl3_route_metadata = _ecp_link_step_route_metadata(sl3_route_certificates)
        return (; from_certificate, to_certificate, endpoint_transport_matrix, factors,
            sl3_route_matrices, sl3_route_certificates, sl3_route_factor_groups,
            sl3_route_metadata)
    end
    throw(ArgumentError("unsupported ECP link step family $(support_family)"))
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
    return _ecp_link_step_endpoint_transport(R, from_column, to_column, link_identity, support_family)
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
    isempty(positions) && return identity
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
    sl3_route_matrices,
    sl3_route_certificates,
    sl3_route_factor_groups,
    sl3_route_metadata,
)
    n = length(from_column)
    from_matrix = matrix(R, n, 1, collect(from_column))
    to_matrix = matrix(R, n, 1, collect(to_column))
    sl2_block_ok = nrows(sl2_block) == 2 &&
        ncols(sl2_block) == 2 &&
        _same_base_ring(base_ring(sl2_block), R) &&
        det(sl2_block) == one(R)
    sl2_embedding_ok = sl2_embedding == block_embedding(sl2_block, n, (1, 2))
    support_family_ok = support_family in (
        :supplied_fixture_identity_sl2_endpoint_transport,
        :direct_elementary_endpoint_transport,
        :polynomial_sl3_route_endpoint_transport,
    )
    endpoint_reductions_ok = verify_ecp_column_reduction(from_certificate) &&
        verify_ecp_column_reduction(to_certificate) &&
        from_certificate.original_column == collect(from_column) &&
        to_certificate.original_column == collect(to_column)
    endpoint_transport_ok = endpoint_transport_matrix == _factor_sequence_product(forward_factors, R, n)
    elementary_factors_ok = _ecp_factor_sequences_equal(elementary_factors, forward_factors)
    route_lengths_ok = length(sl3_route_matrices) == length(sl3_route_certificates) &&
        length(sl3_route_matrices) == length(sl3_route_factor_groups) &&
        length(sl3_route_matrices) == length(sl3_route_metadata)
    route_certificates_ok = route_lengths_ok && all((
        _verify_polynomial_factorization_route_certificate(sl3_route_certificates[idx]) &&
        sl3_route_certificates[idx].status == :supported &&
        sl3_route_certificates[idx].matrix == sl3_route_matrices[idx] &&
        tuple(sl3_route_certificates[idx].factors...) == sl3_route_factor_groups[idx]
        for idx in eachindex(sl3_route_certificates)
    ))
    route_metadata_ok = route_lengths_ok && all((
        get(sl3_route_metadata[idx], :source, nothing) == :polynomial_factorization_route_certificate &&
        get(sl3_route_metadata[idx], :obligation_index, nothing) == idx &&
        get(sl3_route_metadata[idx], :route, nothing) == sl3_route_certificates[idx].route &&
        get(sl3_route_metadata[idx], :status, nothing) == sl3_route_certificates[idx].status &&
        get(sl3_route_metadata[idx], :factor_count, nothing) == length(sl3_route_certificates[idx].factors)
        for idx in eachindex(sl3_route_metadata)
    ))
    route_family_ok = if support_family == :supplied_fixture_identity_sl2_endpoint_transport
        isempty(sl3_route_matrices) &&
            isempty(sl3_route_certificates) &&
            isempty(sl3_route_factor_groups) &&
            isempty(sl3_route_metadata) &&
            sl2_block == identity_matrix(R, 2)
    elseif support_family == :direct_elementary_endpoint_transport
        isempty(sl3_route_matrices) &&
            isempty(sl3_route_certificates) &&
            isempty(sl3_route_factor_groups) &&
            isempty(sl3_route_metadata) &&
            all(factor -> _ecp_is_elementary_factor(factor, R, n), forward_factors)
    else
        route_lengths_ok && route_certificates_ok && route_metadata_ok
    end
    forward_map_ok = _apply_reduction_factors(forward_factors, collect(from_column), R) == to_matrix
    inverse_map_ok = _apply_reduction_factors(inverse_factors, collect(to_column), R) == from_matrix
    inverse_sequence_ok = _ecp_factor_sequences_equal(inverse_factors, _ecp_inverse_factor_sequence(forward_factors))
    overall_ok = sl2_block_ok &&
        sl2_embedding_ok &&
        support_family_ok &&
        endpoint_reductions_ok &&
        endpoint_transport_ok &&
        elementary_factors_ok &&
        route_family_ok &&
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
        route_family_ok,
        route_lengths_ok,
        route_certificates_ok,
        route_metadata_ok,
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
    route_mode_ok = certificate.route_mode ==
        _ecp_link_step_resolve_route_mode(certificate.link_witness, certificate.route_mode)
    transformed_column_ok = !isempty(recomputed_path_columns) &&
        certificate.transformed_column == certificate.original_column &&
        certificate.transformed_column == recomputed_path_columns[end]
    lower_variable_column_ok = !isempty(recomputed_path_columns) &&
        certificate.lower_variable_column == recomputed_path_columns[1]

    recomputed_segments = try
        path_columns_ok ? _ecp_link_step_segments(
            certificate.link_witness,
            recomputed_path_columns;
            route_mode = certificate.route_mode,
        ) : ()
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
        route_mode_ok &&
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
        route_mode_ok,
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
        tuple(left.sl3_route_matrices...) == tuple(right.sl3_route_matrices...) &&
        _ecp_route_certificates_equivalent(left.sl3_route_certificates, right.sl3_route_certificates) &&
        tuple(left.sl3_route_factor_groups...) == tuple(right.sl3_route_factor_groups...) &&
        left.sl3_route_metadata == right.sl3_route_metadata &&
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

function _ecp_route_certificates_equivalent(left, right)
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _verify_polynomial_factorization_route_certificate(left[idx]) || return false
        _verify_polynomial_factorization_route_certificate(right[idx]) || return false
        left[idx].matrix == right[idx].matrix || return false
        left[idx].route == right[idx].route || return false
        tuple(left[idx].factors...) == tuple(right[idx].factors...) || return false
        left[idx].product == right[idx].product || return false
        left[idx].status == right[idx].status || return false
        left[idx].verification == right[idx].verification || return false
    end
    return true
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
            :inverse_substituted_coordinate_move_factors,
            :inverse_substituted_reduction_factors,
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
        source_order_index = findfirst(variable -> variable == source_variable, variable_order)
        target_order_index = findfirst(variable -> variable == target_variable, variable_order)
        summary = _ecp_monicity_normalization_summary(
            collect(input_column),
            R,
            normalized_variable_order,
            target_variable_index;
            source_variable_index,
            shift_power = stage.shift_power,
            shift_sign = stage.shift_sign,
            selected_monic_index_hint = stage.selected_monic_index,
        )
        summary === nothing && return invalid_replay
        transformed_replay = _ecp_replay_stage(
            stage.transformed_stage,
            collect(summary.normalized_column),
            R,
        )
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
            stage.shift_polynomial == summary.shift_polynomial &&
            stage.forward_values == summary.forward_values &&
            stage.inverse_values == summary.inverse_values &&
            stage.forward_substitution == summary.forward_substitution &&
            stage.inverse_substitution == summary.inverse_substitution &&
            stage.transformed_column == summary.transformed_column &&
            stage.selected_monic_index == summary.selected_monic_index &&
            stage.selected_monic_entry == summary.selected_monic_entry &&
            stage.first_coordinate_strategy == summary.first_coordinate_strategy &&
            _ecp_factor_sequences_equal(
                stage.first_coordinate_move_factors,
                summary.coordinate_move_factors,
            ) &&
            stage.first_coordinate_column == summary.normalized_column &&
            transformed_replay.ok &&
            _ecp_factor_sequences_equal(
                stage.transformed_factors,
                summary.transformed_factors,
            ) &&
            _ecp_factor_sequences_equal(
                transformed_replay.factors,
                summary.transformed_reduction_factors,
            ) &&
            _ecp_factor_sequences_equal(
                stage.inverse_substituted_coordinate_move_factors,
                summary.inverse_substituted_coordinate_move_factors,
            ) &&
            _ecp_factor_sequences_equal(
                stage.inverse_substituted_reduction_factors,
                summary.inverse_substituted_reduction_factors,
            ) &&
            _ecp_factor_sequences_equal(
                stage.inverse_substituted_factors,
                summary.inverse_substituted_factors,
            ) &&
            _ecp_factor_sequences_equal(stage.factors, summary.inverse_substituted_factors) &&
            stage.variable_change_verification == summary.variable_change_verification &&
            summary.variable_change_verification.selected_monic_ok &&
            summary.variable_change_verification.substitution_inverse_ok &&
            summary.variable_change_verification.transformed_reduction_ok &&
            summary.variable_change_verification.original_reduction_ok &&
            stage.output_column == summary.original_output
        return (; ok, factors = summary.inverse_substituted_factors, output_column = summary.original_output)
    elseif stage.kind == :ecp_pipeline
        invalid_replay = (; ok = false, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
        required_keys = (
            :kind,
            :input_column,
            :route_metadata,
            :context,
            :normalization,
            :link_witness,
            :link_step,
            :induction_normality,
            :inverse_substituted_induction_factors,
            :factors,
            :output_column,
        )
        _ecp_stage_keys_ok(stage, required_keys) || return invalid_replay

        context_ok = verify_ecp_input_context(stage.context) &&
            tuple(stage.context.column...) == tuple(input_column...) &&
            _same_base_ring(stage.context.ring, R)
        normalization_ok = verify_ecp_monicity_normalization(stage.normalization) &&
            stage.normalization.context == stage.context
        link_witness_ok = verify_ecp_link_witness(stage.link_witness) &&
            stage.link_witness == stage.link_step.link_witness
        link_step_ok = verify_ecp_link_step_certificate(stage.link_step) &&
            stage.link_step.original_column == stage.normalization.normalized_column &&
            _same_base_ring(stage.link_step.ring, R)
        induction_ok = verify_ecp_induction_normality_certificate(stage.induction_normality) &&
            stage.induction_normality.link_step == stage.link_step
        context_ok && normalization_ok && link_witness_ok && link_step_ok && induction_ok || return invalid_replay

        inverse_substituted_induction_factors = _ecp_substitute_factor_sequence(
            stage.induction_normality.final_factors,
            stage.normalization.inverse_substitution,
            R,
        )
        expected_factors = vcat(
            inverse_substituted_induction_factors,
            stage.normalization.inverse_substituted_coordinate_move_factors,
        )
        expected_output = _apply_reduction_factors(expected_factors, input_column, R)
        expected_route_metadata = _ecp_general_pipeline_route_metadata(
            stage.link_step,
            length(stage.normalization.normalized_column),
        )
        ok = stage.input_column == _ecp_column_tuple(input_column) &&
            stage.route_metadata == expected_route_metadata &&
            _ecp_factor_sequences_equal(
                stage.inverse_substituted_induction_factors,
                inverse_substituted_induction_factors,
            ) &&
            _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
            stage.output_column == expected_output
        return (; ok, factors = expected_factors, output_column = expected_output)
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
    elseif stage.kind == :laurent_elementary_row_preconditioning
        invalid_replay = (; ok = false, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
        _ecp_stage_keys_ok(
            stage,
            (
                :kind,
                :input_column,
                :target_index,
                :source_index,
                :source_indices,
                :coefficient,
                :coefficients,
                :coefficient_strategy,
                :precondition_factors,
                :precondition_factor,
                :transformed_column,
                :transformed_certificate,
                :factors,
                :output_column,
            ),
        ) || return invalid_replay
        _is_laurent_polynomial_ring(R) || return invalid_replay

        n = length(input_column)
        target_idx = stage.target_index
        target_idx isa Integer || return invalid_replay
        1 <= target_idx <= n || return invalid_replay

        stage.source_index isa Integer || return invalid_replay
        stage.source_indices isa Tuple || return invalid_replay
        all(source_idx -> source_idx isa Integer, stage.source_indices) ||
            return invalid_replay
        source_indices = stage.source_indices
        isempty(source_indices) && return invalid_replay
        coefficients = tuple((
            _coerce_into_ring(R, coeff, "row preconditioning coefficient")
            for coeff in stage.coefficients
        )...)
        length(source_indices) == length(coefficients) || return invalid_replay
        stage.source_index == source_indices[1] || return invalid_replay
        stage.coefficient == coefficients[1] || return invalid_replay
        all(source_idx -> 1 <= source_idx <= n && source_idx != target_idx, source_indices) ||
            return invalid_replay
        length(unique(source_indices)) == length(source_indices) || return invalid_replay
        _laurent_row_preconditioning_stage_spec_ok(
            input_column,
            R,
            target_idx,
            source_indices,
            coefficients,
            stage.coefficient_strategy,
        ) || return invalid_replay

        precondition_factors = [
            elementary_matrix(n, target_idx, source_idx, coeff, R)
            for (source_idx, coeff) in zip(source_indices, coefficients)
        ]
        precondition_factor = _factor_sequence_product(precondition_factors, R, n)
        transformed_matrix = precondition_factor * matrix(R, n, 1, collect(input_column))
        transformed_column = collect(_ecp_matrix_column_to_tuple(transformed_matrix))
        transformed_certificate_ok =
            verify_ecp_column_reduction(stage.transformed_certificate) &&
            stage.transformed_certificate.original_column == transformed_column &&
            stage.transformed_certificate.ring == R
        expected_factors = transformed_certificate_ok ?
            vcat(stage.transformed_certificate.factors, precondition_factors) :
            Any[]
        expected_output = transformed_certificate_ok ?
            _apply_reduction_factors(expected_factors, input_column, R) :
            matrix(R, n, 1, collect(input_column))
        ok = stage.input_column == _ecp_column_tuple(input_column) &&
            _ecp_factor_sequences_equal(
                stage.precondition_factors,
                precondition_factors,
            ) &&
            stage.precondition_factor == precondition_factor &&
            stage.transformed_column == tuple(transformed_column...) &&
            transformed_certificate_ok &&
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

function _ecp_substitution_map_values(substitution_map)
    return tuple((entry.value for entry in substitution_map)...)
end

function _ecp_substitution_maps_are_inverse(R, forward_substitution, inverse_substitution)::Bool
    try
        ring_gens = tuple(gens(R)...)
        forward_values = _ecp_substitution_map_values(forward_substitution)
        inverse_values = _ecp_substitution_map_values(inverse_substitution)
        length(forward_values) == length(ring_gens) || return false
        length(inverse_values) == length(ring_gens) || return false
        forward_then_inverse = tuple((
            _coerce_into_ring(
                R,
                evaluate(evaluate(generator, collect(forward_values)), collect(inverse_values)),
                "forward-inverse substitution generator",
            )
            for generator in ring_gens
        )...)
        inverse_then_forward = tuple((
            _coerce_into_ring(
                R,
                evaluate(evaluate(generator, collect(inverse_values)), collect(forward_values)),
                "inverse-forward substitution generator",
            )
            for generator in ring_gens
        )...)
        return forward_then_inverse == ring_gens && inverse_then_forward == ring_gens
    catch err
        err isa InterruptException && rethrow()
        return false
    end
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

function _diagnose_laurent_unimodular_column_reduction(
    column::AbstractVector,
    R,
    attempted::Vector{Symbol},
    details::Vector;
    laurent_large_support_diagnostic_decline::Bool = false,
)
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    if unit_idx !== nothing
        push!(details, _column_reduction_stage_detail(:unit_entry, R, :supported; pivot_index = unit_idx))
        return (; supported = true, stage = :unit_entry)
    end
    push!(details, _column_reduction_stage_detail(:unit_entry, R, :no_unit_entry; pivot_index = nothing))

    push!(attempted, :laurent_unit_creation)
    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    if unit_creation !== nothing
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_unit_creation,
                R,
                :supported;
                pivot_index = unit_creation.stage.pivot_index,
                source_index = unit_creation.stage.source_index,
            ),
        )
        return (; supported = true, stage = :laurent_unit_creation)
    end
    push!(details, _column_reduction_stage_detail(:laurent_unit_creation, R, :no_unit_creation_candidate))

    push!(attempted, :laurent_witness_unit)
    witness_decline = laurent_large_support_diagnostic_decline ?
        _laurent_diagnostic_large_support_decline(column) :
        nothing
    witness = if witness_decline === nothing
        _laurent_unimodular_witness(column, R)
    else
        nothing
    end
    if witness === nothing
        outcome = witness_decline === nothing ?
            :witness_unavailable :
            :witness_support_too_large
        witness_decline_detail = witness_decline === nothing ? (;) : witness_decline
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_witness_unit,
                R,
                outcome;
                witness_unit_index = nothing,
                witness_decline_detail...,
            ),
        )
    else
        witness_unit_idx = findfirst(is_unit, witness)
        if witness_unit_idx === nothing
            push!(
                details,
                _column_reduction_stage_detail(
                    :laurent_witness_unit,
                    R,
                    :witness_without_unit;
                    witness_unit_index = nothing,
                ),
            )
        else
            _witness_unit_reduction_certificate_stage(column, witness, witness_unit_idx, R)
            push!(
                details,
                _column_reduction_stage_detail(
                    :laurent_witness_unit,
                    R,
                    :supported;
                    witness_unit_index = witness_unit_idx,
                ),
            )
            return (; supported = true, stage = :laurent_witness_unit)
        end
    end

    push!(attempted, :laurent_normalization)
    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring

    normalized_decline = laurent_large_support_diagnostic_decline ?
        _laurent_diagnostic_large_support_decline(poly_column) :
        nothing
    if normalized_decline !== nothing
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_normalization,
                R,
                :delegation_declined_large_support;
                normalized_column_length = length(poly_column),
                normalized_ring_kind = _column_reduction_ring_kind(P),
                normalized_status = :declined,
                normalized_failure_code = :support_too_large,
                normalized_decline...,
            ),
        )
        return _diagnose_laurent_row_preconditioning(column, R, attempted, details)
    end

    normalized_unimodular = try
        is_unimodular_column(poly_column, P)
    catch err
        err isa InterruptException && rethrow()
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_normalization,
                R,
                :normalized_unimodularity_check_failed;
                normalized_column_length = length(poly_column),
                normalized_ring_kind = _column_reduction_ring_kind(P),
                normalized_status = :precondition_failed,
                normalized_failure_code = :unimodularity_check_failed,
                normalized_message = _column_reduction_error_message(err),
            ),
        )
        return (; supported = false, stage = nothing)
    end

    if !normalized_unimodular
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_normalization,
                R,
                :normalized_not_unimodular;
                normalized_column_length = length(poly_column),
                normalized_ring_kind = _column_reduction_ring_kind(P),
                normalized_status = :precondition_failed,
                normalized_failure_code = :not_unimodular,
            ),
        )
        return _diagnose_laurent_row_preconditioning(column, R, attempted, details)
    end

    ordinary_attempted = Symbol[]
    ordinary_details = Any[]
    result = _diagnose_polynomial_unimodular_column_reduction(
        poly_column,
        P,
        ordinary_attempted,
        ordinary_details;
        allow_general_ecp_pipeline = false,
    )
    normalized_status = result.supported ? :supported : :unsupported
    normalized_failure_code = result.supported ? nothing : :unsupported_polynomial_column_family
    push!(
        details,
        _column_reduction_stage_detail(
            :laurent_normalization,
            R,
            :delegated_to_polynomial;
            normalized_column_length = length(poly_column),
            normalized_ring_kind = _column_reduction_ring_kind(P),
            normalized_status,
            normalized_failure_code,
        ),
    )
    append!(attempted, ordinary_attempted)
    append!(details, ordinary_details)
    result.supported && return result

    return _diagnose_laurent_row_preconditioning(column, R, attempted, details)
end

function _diagnose_laurent_row_preconditioning(
    column::AbstractVector,
    R,
    attempted::Vector{Symbol},
    details::Vector,
)
    push!(attempted, :laurent_elementary_row_preconditioning)
    row_preconditioning = _laurent_row_preconditioning_candidate(column, R)
    if row_preconditioning !== nothing
        transformed_stage = row_preconditioning.transformed_certificate.stages[end].kind
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_elementary_row_preconditioning,
                R,
                :supported;
                target_index = row_preconditioning.target_index,
                source_index = row_preconditioning.source_index,
                source_indices = row_preconditioning.source_indices,
                coefficient = row_preconditioning.coefficient,
                coefficients = row_preconditioning.coefficients,
                coefficient_strategy = row_preconditioning.coefficient_strategy,
                coefficient_count = length(row_preconditioning.coefficients),
                transformed_stage,
            ),
        )
        return (; supported = true, stage = :laurent_elementary_row_preconditioning)
    end
    push!(
        details,
        _column_reduction_stage_detail(
            :laurent_elementary_row_preconditioning,
            R,
            :no_row_preconditioning_candidate;
            target_index = nothing,
            source_index = nothing,
            source_indices = (),
            coefficient = nothing,
            coefficients = (),
            coefficient_strategy = nothing,
            coefficient_count = 0,
            transformed_stage = nothing,
        ),
    )
    push!(attempted, :laurent_native_ecp_boundary)
    push!(details, _laurent_native_ecp_boundary_stage_detail(R))
    return (; supported = false, stage = :laurent_native_ecp_boundary)
end

function _diagnose_polynomial_unimodular_column_reduction(
    column::AbstractVector,
    R,
    attempted::Vector{Symbol},
    details::Vector,
    ;
    allow_general_ecp_pipeline::Bool = true,
)
    supported = _diagnose_supported_unimodular_column_reduction(column, R, attempted, details)
    supported.supported && return supported

    if length(column) > 3
        push!(attempted, :three_entry_block)
        block = _reduce_via_supported_three_block_certificate(column, R)
        if block !== nothing
            push!(
                details,
                _column_reduction_stage_detail(
                    :three_entry_block,
                    R,
                    :supported;
                    block_indices = block.stage.indices,
                    pivot_index = block.stage.indices[end],
                ),
            )
            return (; supported = true, stage = :three_entry_block)
        end
        push!(details, _column_reduction_stage_detail(:three_entry_block, R, :no_supported_three_block))

        if allow_general_ecp_pipeline
            push!(attempted, :general_ecp_pipeline)
            general_failed = false
            general = try
                _reduce_via_general_ecp_pipeline_certificate(column, R)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
                general_failed = true
                push!(
                    details,
                    _column_reduction_stage_detail(
                        :general_ecp_pipeline,
                        R,
                        :staged_failure;
                        message = _column_reduction_error_message(err),
                    ),
                )
                nothing
            end
            if general !== nothing
                push!(
                    details,
                    _column_reduction_stage_detail(
                        :general_ecp_pipeline,
                        R,
                        :supported;
                        link_route_mode = general.stage.route_metadata.link_route_mode,
                        normalized_column_length = general.stage.route_metadata.normalized_column_length,
                    ),
                )
                return (; supported = true, stage = :general_ecp_pipeline)
            end
            general_failed ||
                push!(details, _column_reduction_stage_detail(:general_ecp_pipeline, R, :unsupported))
        end

        push!(attempted, :monicity_normalization)
        normalized = try
            _reduce_after_monicity_normalization_certificate(column, R)
        catch err
            err isa InterruptException && rethrow() # COV_EXCL_LINE
            nothing
        end
        if normalized !== nothing
            push!(
                details,
                _column_reduction_stage_detail(
                    :monicity_normalization,
                    R,
                    :supported;
                    normalized_column_length = length(column),
                ),
            )
            return (; supported = true, stage = :monicity_normalization)
        end
        push!(details, _column_reduction_stage_detail(:monicity_normalization, R, :no_monicity_normalization))
        return (; supported = false, stage = nothing)
    end

    if allow_general_ecp_pipeline
        push!(attempted, :general_ecp_pipeline)
        general_failed = false
        general = try
            _reduce_via_general_ecp_pipeline_certificate(column, R)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            general_failed = true
            push!(
                details,
                _column_reduction_stage_detail(
                    :general_ecp_pipeline,
                    R,
                    :staged_failure;
                    message = _column_reduction_error_message(err),
                ),
            )
            nothing
        end
        if general !== nothing
            push!(
                details,
                _column_reduction_stage_detail(
                    :general_ecp_pipeline,
                    R,
                    :supported;
                    link_route_mode = general.stage.route_metadata.link_route_mode,
                    normalized_column_length = general.stage.route_metadata.normalized_column_length,
                ),
            )
            return (; supported = true, stage = :general_ecp_pipeline)
        end
        general_failed ||
            push!(details, _column_reduction_stage_detail(:general_ecp_pipeline, R, :unsupported))
    end

    push!(attempted, :monicity_normalization)
    normalized = try
        _reduce_after_monicity_normalization_certificate(column, R)
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        nothing
    end
    if normalized !== nothing
        push!(
            details,
            _column_reduction_stage_detail(
                :monicity_normalization,
                R,
                :supported;
                normalized_column_length = length(column),
            ),
        )
        return (; supported = true, stage = :monicity_normalization)
    end
    push!(details, _column_reduction_stage_detail(:monicity_normalization, R, :no_monicity_normalization))

    return (; supported = false, stage = nothing)
end

function _diagnose_exact_small_column_reduction(
    column::AbstractVector,
    R,
    attempted::Vector{Symbol},
    details::Vector,
)
    supported = _diagnose_supported_unimodular_column_reduction(column, R, attempted, details)
    supported.supported && return supported

    if _has_at_least_two_generators(R)
        push!(attempted, :monicity_normalization)
        normalized = try
            _reduce_after_monicity_normalization_certificate(column, R)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        if normalized !== nothing
            push!(
                details,
                _column_reduction_stage_detail(
                    :monicity_normalization,
                    R,
                    :supported;
                    normalized_column_length = length(column),
                ),
            )
            return (; supported = true, stage = :monicity_normalization)
        end
        push!(details, _column_reduction_stage_detail(:monicity_normalization, R, :no_monicity_normalization))
    end

    return (; supported = false, stage = nothing)
end

function _diagnose_supported_unimodular_column_reduction(
    column::AbstractVector,
    R,
    attempted::Vector{Symbol},
    details::Vector,
)
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    if unit_idx !== nothing
        push!(details, _column_reduction_stage_detail(:unit_entry, R, :supported; pivot_index = unit_idx))
        return (; supported = true, stage = :unit_entry)
    end
    push!(details, _column_reduction_stage_detail(:unit_entry, R, :no_unit_entry; pivot_index = nothing))

    push!(attempted, :witness_unit)
    witness = try
        _unimodular_witness(column, R)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    if witness === nothing
        push!(
            details,
            _column_reduction_stage_detail(
                :witness_unit,
                R,
                :witness_unavailable;
                witness_unit_index = nothing,
            ),
        )
        return (; supported = false, stage = nothing)
    end
    witness_unit_idx = findfirst(is_unit, witness)
    if witness_unit_idx !== nothing
        push!(
            details,
            _column_reduction_stage_detail(
                :witness_unit,
                R,
                :supported;
                witness_unit_index = witness_unit_idx,
            ),
        )
        return (; supported = true, stage = :witness_unit)
    end

    push!(
        details,
        _column_reduction_stage_detail(
            :witness_unit,
            R,
            :witness_without_unit;
            witness_unit_index = nothing,
        ),
    )

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

function _ecp_first_coordinate_move_factors(n::Int, selected_monic_index::Int, R)
    selected_monic_index == 1 && return typeof(identity_matrix(R, n))[]

    return [
        elementary_matrix(n, 1, selected_monic_index, one(R), R),
        elementary_matrix(n, selected_monic_index, 1, -one(R), R),
        elementary_matrix(n, 1, selected_monic_index, one(R), R),
    ]
end

function _ecp_first_coordinate_strategy(selected_monic_index::Int)
    return selected_monic_index == 1 ? :already_first : :moved_to_first
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

function _unsupported_unimodular_column_reduction_message(
    column::AbstractVector,
    R;
    detail = nothing,
)
    profile = _is_laurent_polynomial_ring(R) ? "Laurent-normalized" : "ordinary polynomial"
    n = length(column)
    prefix = "unsupported exact unimodular column reduction for $(profile) column of length $(n): "
    if detail !== nothing
        return prefix * detail
    end
    if _is_laurent_polynomial_ring(R)
        return prefix *
            "no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies"
    end
    return prefix *
        "no supported unit, witness-unit, monicity-normalized, 3-entry block, or general ECP pipeline stage applies"
end

function _throw_unsupported_unimodular_column_reduction(column::AbstractVector, R)
    throw(ArgumentError(_unsupported_unimodular_column_reduction_message(column, R)))
end
