struct PolynomialFactorizationRouteCertificate
    matrix
    route::Symbol
    factors::Vector
    product
    evidence
    status::Symbol
    verification
end

struct SL3RealizationInputContext
    matrix
    base_ring
    coefficient_ring
    size::Int
    ring_profile::Symbol
    generators::Tuple
    generator_names::Tuple
    selected_variable
    selected_variable_index
    selected_variable_status::Symbol
    determinant
    determinant_status::Symbol
    exact_field_status::Symbol
    catalog_metadata::NamedTuple
    local_form_witness
    local_form_status::Symbol
    variable_change_metadata
    variable_change_status::Symbol
    normality_conjugation_metadata
    normality_conjugation_status::Symbol
    quillen_murthy_metadata
    quillen_murthy_status::Symbol
    evidence_status::Symbol
    support_status::Symbol
    staged_diagnostic::NamedTuple
    verification
end

struct PolynomialQuillenPatchRouteAdapter
    target
    route::Symbol
    quillen_patch
    global_elementary_factors::Vector
    product
    target_matrix
    replay_metadata
    verification
end

const _POLYNOMIAL_FACTORIZATION_ROUTE_TAGS = Set([
    :fast_local_sl3,
    :disjoint_local_blocks,
    :quillen_patch,
    :polynomial_column_peel,
    :recursive_column_peel,
    :staged_failure,
])

_is_polynomial_column_peel_route(route::Symbol) =
    route in (:polynomial_column_peel, :recursive_column_peel)

function _validate_factorization_matrix(A)
    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))
    n = nrows(A)
    n >= 3 || throw(ArgumentError("elementary_factorization requires matrices of size at least 3"))
    return n
end

function _factorization_ring_profile(R)
    _is_laurent_polynomial_ring(R) && return :laurent
    try
        collect(gens(R))
        return :polynomial
    catch err
        err isa MethodError || rethrow()
        throw(ArgumentError("A base ring is outside the supported exact polynomial or Laurent polynomial factorization path"))
    end
end

function _polynomial_exact_field_backed_ring(R)::Bool
    try
        return Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) &&
               coefficient_ring(R) isa Field
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_unsupported_coefficient_ring_message()
    return "ordinary polynomial factorization currently requires an exact " *
           "field-backed coefficient ring; coefficient-ring support beyond " *
           "exact field-backed ordinary polynomial rings remains staged"
end

function _normalize_factorization_input(A, ring_profile::Symbol)
    if ring_profile == :laurent
        normalization = normalize_laurent_gl_matrix(A)
        return normalization.normalized_matrix, normalization
    end

    return A, nothing
end

function _require_polynomial_sl_determinant(A)
    R = base_ring(A)
    det(A) == one(R) || throw(ArgumentError("determinant/unit precondition failed: polynomial inputs must have determinant 1; otherwise the input is outside the staged SL_n factorization path"))
    return nothing
end

function _sl3_realization_input_context_extract(data, fields::Tuple)
    data === nothing && return nothing
    for field in fields
        hasproperty(data, field) && return getproperty(data, field)
    end
    return nothing
end

function _sl3_realization_input_context_nonempty_identifier(value)::Bool
    value === nothing && return false
    value isa AbstractString && return !isempty(value)
    return true
end

function _sl3_realization_input_context_has_replay_payload(metadata)::Bool
    payload = _sl3_realization_input_context_extract(
        metadata,
        (
            :replay_steps,
            :replay_metadata,
            :replay_certificate,
            :replay_payload,
            :certificate,
            :verification,
            :patch,
            :variable_change_certificate,
            :normality_certificate,
            :conjugation_certificate,
            :quillen_patch,
            :murthy_certificate,
        ),
    )
    payload === nothing && return false
    if payload isa Tuple || payload isa AbstractArray || payload isa Set
        return !isempty(payload)
    end
    return true
end

function _sl3_realization_input_context_selected_variable(R, selected_variable)
    selected_variable === nothing && return nothing, nothing, :missing

    normalized_selected =
        hasproperty(selected_variable, :generator) ? getproperty(selected_variable, :generator) :
        selected_variable
    index = _ecp_selected_variable_index(R, normalized_selected)

    if hasproperty(selected_variable, :index)
        metadata_index = getproperty(selected_variable, :index)
        metadata_index == index || throw(ArgumentError(
            "selected_variable metadata index does not match the generator position in gens(R)",
        ))
    end

    return _require_substitution_generator(R, normalized_selected), index, :passes
end

function _sl3_realization_input_context_local_form_status(
    A,
    R,
    ring_profile::Symbol,
    selected_variable,
    selected_variable_index,
    local_form_witness,
)
    supported_generator = _supported_local_sl3_generator(A, R, ring_profile)
    if local_form_witness !== nothing
        monic_value = _sl3_realization_input_context_extract(
            local_form_witness,
            (:entry, :polynomial, :p, :local_entry),
        )
        position = hasproperty(local_form_witness, :monic_entry_position) ?
            local_form_witness.monic_entry_position :
            nothing
        if monic_value !== nothing &&
                position isa Tuple &&
                length(position) == 2 &&
                selected_variable_index !== nothing
            row, col = position
            if row isa Int &&
                    col isa Int &&
                    1 <= row <= nrows(A) &&
                    1 <= col <= ncols(A) &&
                    parent(monic_value) === R &&
                    A[row, col] == monic_value &&
                    _is_monic_in_variable(monic_value, R, selected_variable_index)
                return :replayed
            end
        end
    end

    if local_form_witness === nothing &&
            supported_generator !== nothing &&
            selected_variable == supported_generator
        return :fast_local
    end

    return :missing
end

function _sl3_realization_input_context_metadata_status(metadata, A)
    metadata === nothing && return :missing

    replay_id = _sl3_realization_input_context_extract(
        metadata,
        (:replay_id, :case_id, :mainline_case_id, :upstream_case_id, :patch_case_id, :fixture_id),
    )
    target_matrix = _sl3_realization_input_context_extract(
        metadata,
        (:target_matrix, :expected_matrix, :matrix),
    )
    has_replay_payload = _sl3_realization_input_context_has_replay_payload(metadata)
    has_replay_id = _sl3_realization_input_context_nonempty_identifier(replay_id)

    if has_replay_id && target_matrix !== nothing && target_matrix == A && has_replay_payload
        return :replayed
    elseif has_replay_id || target_matrix !== nothing || has_replay_payload
        return :recorded
    end

    return :recorded
end

function _sl3_realization_input_context_staged_diagnostic(
    selected_variable_status::Symbol,
    determinant_status::Symbol,
    exact_field_status::Symbol,
    local_form_status::Symbol,
    variable_change_status::Symbol,
    normality_conjugation_status::Symbol,
    quillen_murthy_status::Symbol,
    evidence_status::Symbol,
    support_status::Symbol,
)
    missing_evidence = Symbol[]
    partial_evidence = Symbol[]

    local_form_status in (:replayed, :fast_local) || push!(missing_evidence, :local_form)
    for (label, status) in (
        (:variable_change, variable_change_status),
        (:normality_conjugation, normality_conjugation_status),
        (:quillen_murthy, quillen_murthy_status),
    )
        if status == :recorded
            push!(partial_evidence, label)
        elseif status != :replayed
            push!(missing_evidence, label)
        end
    end

    message = support_status == :supported ?
        "SL_3 realization input context is replayable" :
        evidence_status == :partial ?
        "SL_3 realization input context is staged with partial replay metadata" :
        "SL_3 realization input context is staged: replayable evidence is missing"

    return (;
        status = support_status,
        missing_evidence = Tuple(missing_evidence),
        partial_evidence = Tuple(partial_evidence),
        selected_variable_status,
        determinant_status,
        exact_field_status,
        message,
    )
end

function _sl3_realization_input_context_fields(
    A;
    selected_variable = nothing,
    catalog_metadata = (;),
    local_form_witness = nothing,
    variable_change_metadata = nothing,
    normality_conjugation_metadata = nothing,
    quillen_murthy_metadata = nothing,
)
    size = _validate_factorization_matrix(A)
    size == 3 || throw(ArgumentError("SL_3 realization input context requires a 3 x 3 matrix"))

    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    ring_profile == :polynomial ||
        throw(ArgumentError("SL_3 realization input context requires an ordinary polynomial base ring"))

    coefficient = coefficient_ring(R)
    Oscar.is_exact_type(typeof(zero(coefficient))) ||
        throw(ArgumentError("SL_3 realization input context requires an exact coefficient ring"))
    coefficient isa Field ||
        throw(ArgumentError("SL_3 realization input context requires a field-backed coefficient ring"))

    _require_polynomial_sl_determinant(A)
    determinant = det(A)

    generators = Tuple(collect(gens(R)))
    generator_names = Tuple(string(generator) for generator in generators)
    selected, selected_index, selected_status =
        _sl3_realization_input_context_selected_variable(R, selected_variable)

    local_form_status = _sl3_realization_input_context_local_form_status(
        A,
        R,
        ring_profile,
        selected,
        selected_index,
        local_form_witness,
    )
    variable_change_status = _sl3_realization_input_context_metadata_status(variable_change_metadata, A)
    normality_conjugation_status =
        _sl3_realization_input_context_metadata_status(normality_conjugation_metadata, A)
    quillen_murthy_status = _sl3_realization_input_context_metadata_status(quillen_murthy_metadata, A)

    replayable_evidence =
        local_form_status in (:replayed, :fast_local) ||
        variable_change_status == :replayed ||
        normality_conjugation_status == :replayed ||
        quillen_murthy_status == :replayed
    recorded_only_evidence =
        !replayable_evidence &&
        (
            variable_change_status == :recorded ||
            normality_conjugation_status == :recorded ||
            quillen_murthy_status == :recorded
        )

    evidence_status = replayable_evidence ? :replayable : recorded_only_evidence ? :partial : :missing
    support_status = replayable_evidence ? :supported : :staged
    staged_diagnostic = _sl3_realization_input_context_staged_diagnostic(
        selected_status,
        :one,
        :supported,
        local_form_status,
        variable_change_status,
        normality_conjugation_status,
        quillen_murthy_status,
        evidence_status,
        support_status,
    )

    if hasproperty(catalog_metadata, :expected_status)
        catalog_metadata.expected_status == support_status || throw(
            ArgumentError("SL_3 realization input context support status does not match catalog metadata")
        )
    end

    return (;
        matrix = A,
        base_ring = R,
        coefficient_ring = coefficient,
        size,
        ring_profile,
        generators,
        generator_names,
        selected_variable = selected,
        selected_variable_index = selected_index,
        selected_variable_status = selected_status,
        determinant,
        determinant_status = :one,
        exact_field_status = :supported,
        catalog_metadata,
        local_form_witness,
        local_form_status,
        variable_change_metadata,
        variable_change_status,
        normality_conjugation_metadata,
        normality_conjugation_status,
        quillen_murthy_metadata,
        quillen_murthy_status,
        evidence_status,
        support_status,
        staged_diagnostic,
    )
end

function _sl3_realization_input_context_core_verification(context)
    recomputed = _sl3_realization_input_context_fields(
        context.matrix;
        selected_variable = context.selected_variable,
        catalog_metadata = context.catalog_metadata,
        local_form_witness = context.local_form_witness,
        variable_change_metadata = context.variable_change_metadata,
        normality_conjugation_metadata = context.normality_conjugation_metadata,
        quillen_murthy_metadata = context.quillen_murthy_metadata,
    )

    matrix_ok = context.matrix == recomputed.matrix
    base_ring_ok = context.base_ring == recomputed.base_ring
    coefficient_ring_ok = context.coefficient_ring == recomputed.coefficient_ring
    size_ok = context.size == recomputed.size
    ring_profile_ok = context.ring_profile == recomputed.ring_profile
    generators_ok = context.generators == recomputed.generators
    generator_names_ok = context.generator_names == recomputed.generator_names
    selected_variable_ok = context.selected_variable == recomputed.selected_variable
    selected_variable_index_ok =
        context.selected_variable_index == recomputed.selected_variable_index
    selected_variable_status_ok =
        context.selected_variable_status == recomputed.selected_variable_status
    determinant_ok = context.determinant == recomputed.determinant
    determinant_status_ok = context.determinant_status == recomputed.determinant_status
    exact_field_status_ok = context.exact_field_status == recomputed.exact_field_status
    catalog_metadata_ok = context.catalog_metadata == recomputed.catalog_metadata
    local_form_witness_ok = context.local_form_witness == recomputed.local_form_witness
    local_form_status_ok = context.local_form_status == recomputed.local_form_status
    variable_change_metadata_ok =
        context.variable_change_metadata == recomputed.variable_change_metadata
    variable_change_status_ok =
        context.variable_change_status == recomputed.variable_change_status
    normality_conjugation_metadata_ok =
        context.normality_conjugation_metadata == recomputed.normality_conjugation_metadata
    normality_conjugation_status_ok =
        context.normality_conjugation_status == recomputed.normality_conjugation_status
    quillen_murthy_metadata_ok =
        context.quillen_murthy_metadata == recomputed.quillen_murthy_metadata
    quillen_murthy_status_ok =
        context.quillen_murthy_status == recomputed.quillen_murthy_status
    evidence_status_ok = context.evidence_status == recomputed.evidence_status
    support_status_ok = context.support_status == recomputed.support_status
    staged_diagnostic_ok = context.staged_diagnostic == recomputed.staged_diagnostic

    overall_core_ok =
        matrix_ok &&
        base_ring_ok &&
        coefficient_ring_ok &&
        size_ok &&
        ring_profile_ok &&
        generators_ok &&
        generator_names_ok &&
        selected_variable_ok &&
        selected_variable_index_ok &&
        selected_variable_status_ok &&
        determinant_ok &&
        determinant_status_ok &&
        exact_field_status_ok &&
        catalog_metadata_ok &&
        local_form_witness_ok &&
        local_form_status_ok &&
        variable_change_metadata_ok &&
        variable_change_status_ok &&
        normality_conjugation_metadata_ok &&
        normality_conjugation_status_ok &&
        quillen_murthy_metadata_ok &&
        quillen_murthy_status_ok &&
        evidence_status_ok &&
        support_status_ok &&
        staged_diagnostic_ok

    return (;
        matrix_ok,
        base_ring_ok,
        coefficient_ring_ok,
        size_ok,
        ring_profile_ok,
        generators_ok,
        generator_names_ok,
        selected_variable_ok,
        selected_variable_index_ok,
        selected_variable_status_ok,
        determinant_ok,
        determinant_status_ok,
        exact_field_status_ok,
        catalog_metadata_ok,
        local_form_witness_ok,
        local_form_status_ok,
        variable_change_metadata_ok,
        variable_change_status_ok,
        normality_conjugation_metadata_ok,
        normality_conjugation_status_ok,
        quillen_murthy_metadata_ok,
        quillen_murthy_status_ok,
        evidence_status_ok,
        support_status_ok,
        staged_diagnostic_ok,
        overall_core_ok,
    )
end

function _sl3_realization_input_context_verification(context)
    core = _sl3_realization_input_context_core_verification(context)
    stored_verification_ok = context.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_sl3_realization_input_context(context)::Bool
    try
        return _sl3_realization_input_context_verification(context).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sl3_realization_input_context(
    A;
    selected_variable = nothing,
    catalog_metadata = (;),
    local_form_witness = nothing,
    variable_change_metadata = nothing,
    normality_conjugation_metadata = nothing,
    quillen_murthy_metadata = nothing,
)
    fields = _sl3_realization_input_context_fields(
        A;
        selected_variable,
        catalog_metadata,
        local_form_witness,
        variable_change_metadata,
        normality_conjugation_metadata,
        quillen_murthy_metadata,
    )
    context = SL3RealizationInputContext(values(merge(fields, (; verification = nothing,)))...)
    verification = _sl3_realization_input_context_core_verification(context)
    checked = SL3RealizationInputContext(values(merge(fields, (; verification,)))...)
    _verify_sl3_realization_input_context(checked) ||
        error("internal SL_3 realization input context verification failed")
    return checked
end

function _supported_local_sl3_generator(A, R, ring_profile::Symbol)
    ring_profile == :polynomial || return nothing
    nrows(A) == 3 || return nothing

    ring_gens = collect(gens(R))
    length(ring_gens) == 1 || return nothing

    if A[1, 3] != zero(R) || A[2, 3] != zero(R) ||
            A[3, 1] != zero(R) || A[3, 2] != zero(R) ||
            A[3, 3] != one(R)
        return nothing
    end

    return ring_gens[1]
end

function _throw_staged_factorization_failure(A, ring_profile::Symbol, normalization)
    n = nrows(A)

    if ring_profile == :laurent
        classification = normalization.determinant_classification
        throw(ArgumentError("Laurent GL_n normalization boundary succeeded with determinant classification $(classification), but the determinant-correction/driver path cannot yet return elementary factors that reconstruct the original input"))
    end

    if !_polynomial_exact_field_backed_ring(base_ring(A))
        throw(ArgumentError(_polynomial_unsupported_coefficient_ring_message()))
    end

    length(collect(gens(base_ring(A)))) > 1 &&
        throw(ArgumentError("determinant-one polynomial input is outside the implemented Quillen/local evidence route: missing Quillen/local realizability witness"))

    if n > 3
        throw(ArgumentError("SL_n reduction layer for n > 3 is not yet implemented in elementary_factorization"))
    end

    throw(ArgumentError("staged reduction to the supported univariate local SL_3 slice is not yet implemented in elementary_factorization"))
end

function _polynomial_factorization_route_certificate(
    A;
    route=nothing,
    quillen_patch=nothing,
    allow_recursive_column_peel::Bool=true,
)
    n = _validate_factorization_matrix(A)
    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    ring_profile == :polynomial ||
        throw(ArgumentError("polynomial route certificates are only supported for ordinary polynomial inputs"))
    _require_polynomial_sl_determinant(A)
    if !_polynomial_exact_field_backed_ring(R)
        if route === nothing || route == :staged_failure
            return _polynomial_staged_failure_route_certificate(
                A;
                allow_recursive_column_peel = allow_recursive_column_peel,
            )
        end
        throw(ArgumentError(_polynomial_unsupported_coefficient_ring_message()))
    end

    if route === nothing
        quillen_patch === nothing ||
            throw(ArgumentError("Quillen patch route certificates require route = :quillen_patch"))

        X = _supported_local_sl3_generator(A, R, ring_profile)
        X !== nothing && return _polynomial_fast_local_sl3_route_certificate(A, X)

        if n > 3
            try
                return _polynomial_disjoint_local_blocks_route_certificate(A)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
            end
        end

        quillen_certificate = _polynomial_quillen_supplied_evidence_route_certificate(A)
        quillen_certificate !== nothing && return quillen_certificate

        if allow_recursive_column_peel
            try
                return _polynomial_recursive_column_peel_route_certificate(
                    A;
                    route_tag = :polynomial_column_peel,
                )
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
            end
        end

        return _polynomial_staged_failure_route_certificate(
            A;
            allow_recursive_column_peel = allow_recursive_column_peel,
        )
    end

    route isa Symbol || throw(ArgumentError("polynomial route certificate route must be a Symbol"))
    if route == :fast_local_sl3
        X = _supported_local_sl3_generator(A, R, ring_profile)
        X !== nothing ||
            throw(ArgumentError("fast local SL_3 route requires a supported univariate local SL_3 input"))
        return _polynomial_fast_local_sl3_route_certificate(A, X)
    elseif route == :disjoint_local_blocks
        n > 3 ||
            throw(ArgumentError("disjoint local-block route requires a matrix of size greater than 3"))
        return _polynomial_disjoint_local_blocks_route_certificate(A)
    elseif route == :quillen_patch
        quillen_patch !== nothing ||
            throw(ArgumentError("Quillen patch route requires a supplied verified patch"))
        return _polynomial_quillen_patch_route_certificate(A, quillen_patch)
    elseif _is_polynomial_column_peel_route(route)
        return _polynomial_recursive_column_peel_route_certificate(A; route_tag = route)
    elseif route == :staged_failure
        return _polynomial_staged_failure_route_certificate(A)
    end

    throw(ArgumentError("unsupported polynomial factorization route certificate tag $(route)"))
end

function _polynomial_fast_local_sl3_route_certificate(A, X)
    evidence = realize_sl3_local_certificate(A, X)
    factors = copy(evidence.factors)
    product = _polynomial_route_factor_product(factors, base_ring(A), nrows(A))
    return _polynomial_route_certificate(A, :fast_local_sl3, factors, product, evidence, :supported)
end

function _polynomial_disjoint_local_blocks_route_certificate(A)
    evidence = reduce_sln_to_sl3(A)
    factors = copy(evidence.factors)
    product = _polynomial_route_factor_product(factors, base_ring(A), nrows(A))
    return _polynomial_route_certificate(A, :disjoint_local_blocks, factors, product, evidence, :supported)
end

function _polynomial_staged_failure_route_certificate(
    A;
    allow_recursive_column_peel::Bool = true,
)
    R = base_ring(A)
    n = nrows(A)
    product = identity_matrix(R, n)
    factors = typeof(product)[]
    evidence = _polynomial_staged_failure_evidence(
        A;
        allow_recursive_column_peel = allow_recursive_column_peel,
    )
    return _polynomial_route_certificate(A, :staged_failure, factors, product, evidence, :staged)
end

function _polynomial_recursive_column_peel_route_certificate(
    A;
    route_tag::Symbol = :polynomial_column_peel,
)
    _is_polynomial_column_peel_route(route_tag) ||
        throw(ArgumentError("unsupported polynomial column-peel route tag $(route_tag)"))
    evidence = _polynomial_column_peel_certificate(A)
    factors = copy(evidence.factors)
    product = _polynomial_route_factor_product(factors, base_ring(A), nrows(A))
    return _polynomial_route_certificate(A, route_tag, factors, product, evidence, :supported)
end

function _polynomial_quillen_patch_route_certificate(A, patch)
    adapter = _polynomial_quillen_patch_route_adapter(A, patch)
    factors = copy(adapter.global_elementary_factors)
    return _polynomial_route_certificate(
        adapter.target_matrix,
        :quillen_patch,
        factors,
        adapter.product,
        adapter,
        :supported,
    )
end

function _polynomial_quillen_supplied_evidence_route_certificate(A)
    patch = _polynomial_quillen_supplied_evidence_patch(A)
    patch === nothing && return nothing
    return _polynomial_quillen_patch_route_certificate(A, patch)
end

function _polynomial_quillen_elementary_entry(A)
    nrows(A) == 3 && ncols(A) == 3 || return nothing
    R = base_ring(A)
    _factorization_ring_profile(R) == :polynomial || return nothing
    _polynomial_exact_field_backed_ring(R) || return nothing
    ring_gens = collect(gens(R))
    length(ring_gens) >= 2 || return nothing
    row = 0
    col = 0
    entry = zero(R)
    for i in 1:3, j in 1:3
        if i == j
            A[i, j] == one(R) || return nothing
        elseif A[i, j] != zero(R)
            row == 0 || return nothing
            row = i
            col = j
            entry = A[i, j]
        end
    end
    row == 0 && return nothing
    return (;
        row,
        col,
        entry,
        selected_variable = ring_gens[1],
        cover_generator = ring_gens[2],
    )
end

function _polynomial_quillen_supplied_evidence_data(A)
    elementary_entry = _polynomial_quillen_elementary_entry(A)
    elementary_entry === nothing && return nothing

    R = base_ring(A)
    ring_gens = collect(gens(R))
    selected_variable = elementary_entry.selected_variable
    base_entry = evaluate(
        elementary_entry.entry,
        [gen == selected_variable ? zero(R) : gen for gen in ring_gens],
    )
    delta_entry = elementary_entry.entry - base_entry
    delta_entry == zero(R) && return nothing

    factor_type = typeof(identity_matrix(R, 3))
    base_term_policy = base_entry == zero(R) ? :trivial : :supplied
    base_term_factors = base_term_policy == :trivial ?
        factor_type[] :
        [elementary_matrix(3, elementary_entry.row, elementary_entry.col, base_entry, R)]

    denominators = [elementary_entry.cover_generator, one(R) - elementary_entry.cover_generator]
    local_certificates = map(enumerate(denominators)) do (local_index, denominator)
        local_certificate = LocalCertificate(
            [elementary_entry.row, elementary_entry.col],
            [denominator, denominator],
        )
        local_factor = elementary_matrix(
            3,
            elementary_entry.row,
            elementary_entry.col,
            denominator * delta_entry,
            R,
        )
        local_realization = quillen_local_realization_certificate(
            A,
            selected_variable;
            local_certificate = local_certificate,
            denominator = denominator,
            coverage_multiplier = one(R),
            correction = QuillenElementaryCorrection(
                elementary_entry.row,
                elementary_entry.col,
                delta_entry,
            ),
            factors = [local_factor],
            local_correction = local_factor,
            witness_metadata = (;
                source = :automatic_quillen_supplied_evidence,
                local_index = local_index,
                consumer_issue_id = "#220",
            ),
            ring = R,
            size = 3,
        )
        return quillen_local_factor_sequence_certificate(
            local_realization;
            factor_provenance = (;
                factor_index = 1,
                sequence_index = 1,
                local_index = 1,
                source = :automatic_quillen_supplied_evidence,
                consumer_issue_id = "#220",
            ),
            metadata = (;
                source = :automatic_quillen_supplied_evidence,
                local_index = local_index,
                consumer_issue_id = "#220",
            ),
        )
    end

    return (;
        selected_variable,
        cover_generator = elementary_entry.cover_generator,
        row = elementary_entry.row,
        col = elementary_entry.col,
        entry = elementary_entry.entry,
        base_entry,
        delta_entry,
        base_term_policy,
        base_term_factors,
        local_certificates,
    )
end

function _polynomial_quillen_supplied_evidence_patch(A)
    data = _polynomial_quillen_supplied_evidence_data(A)
    data === nothing && return nothing

    R = base_ring(A)
    patch = assemble_quillen_patch_from_local_evidence(
        A,
        data.selected_variable,
        data.local_certificates;
        exponent = 1,
        coverage_multipliers = [one(R), one(R)],
        base_term_policy = data.base_term_policy,
        base_term_factors = data.base_term_factors,
        metadata = (; source = :automatic_quillen_supplied_evidence, consumer_issue_id = "#220"),
    )
    verify_quillen_patch(patch) ||
        error("internal supplied-evidence Quillen patch verification failed")
    return patch
end

function _polynomial_quillen_fixture_data(A)
    nrows(A) == 3 || return nothing
    ncols(A) == 3 || return nothing

    R = base_ring(A)
    _factorization_ring_profile(R) == :polynomial || return nothing
    coefficient_ring(R) == QQ || return nothing

    ring_gens = collect(gens(R))
    length(ring_gens) == 3 || return nothing
    string.(ring_gens) == ["X", "r", "g"] || return nothing

    X, r, g = ring_gens
    entry = X + r^2 * g + g + one(R)
    A == elementary_matrix(3, 1, 2, entry, R) || return nothing

    base_matrix = elementary_matrix(3, 1, 2, X + g + one(R), R)
    witness = (;
        matrix = base_matrix,
        variable = X,
        denominator = r,
        exponent = 2,
        shift = g,
        expected_matrix = A,
    )

    return (;
        ring = R,
        size = 3,
        X,
        r,
        g,
        entry,
        witness,
    )
end

function _polynomial_quillen_fixture_local_certificate(
    A,
    data;
    denominator,
    local_index::Int,
)
    R = data.ring
    n = data.size
    local_factor = elementary_matrix(n, 1, 2, denominator * data.entry, R)
    return quillen_local_realization_certificate(
        A,
        data.X;
        local_certificate = LocalCertificate([1, 2], [denominator, denominator]),
        denominator = denominator,
        coverage_multiplier = one(R),
        correction = QuillenElementaryCorrection(1, 2, data.entry),
        factors = [local_factor],
        local_correction = local_factor,
        patched_substitution_witness = data.witness,
        witness_metadata = (;
            fixture_id = "quillen-patched-substitution-witness-qq",
            local_index = local_index,
            source_refs = ("Issue 99 quillen patched substitution witness shape",),
            consumer_issue_ids = ("#64", "#99", "#105", "#109", "#115"),
        ),
    )
end

_polynomial_quillen_patch_factors(patch) =
    hasproperty(patch, :global_elementary_factors) ? patch.global_elementary_factors :
    hasproperty(patch, :factors) ? patch.factors :
    throw(ArgumentError("Quillen patch has no global elementary factors"))

_polynomial_quillen_patch_product(patch) =
    hasproperty(patch, :patched_product) ? patch.patched_product :
    hasproperty(patch, :product) ? patch.product :
    throw(ArgumentError("Quillen patch has no patched product"))

function _polynomial_quillen_route_target_matrix(target, patch)
    R = _require_supported_quillen_ring(patch.ring)
    _factorization_ring_profile(R) == :polynomial ||
        throw(ArgumentError("Quillen patch route target must lie over an ordinary polynomial ring"))
    n = patch.size
    target_matrix = _quillen_global_target_matrix(
        target;
        ring = R,
        size = n,
        label = "Quillen patch route target",
    )
    _require_square_matrix(target_matrix, "Quillen patch route target") == n ||
        throw(DimensionMismatch("Quillen patch route target size must match the patch size"))
    _same_base_ring(base_ring(target_matrix), R) ||
        throw(ArgumentError("Quillen patch route target must lie in the patch ring"))
    _factorization_ring_profile(base_ring(target_matrix)) == :polynomial ||
        throw(ArgumentError("Quillen patch route target must lie over an ordinary polynomial ring"))
    _require_polynomial_sl_determinant(target_matrix)
    return target_matrix
end

function _polynomial_quillen_patch_target_matrix(patch)
    R = _require_supported_quillen_ring(patch.ring)
    return _quillen_global_target_matrix(
        patch.target;
        ring = R,
        size = patch.size,
        label = "Quillen patch target",
    )
end

function _polynomial_quillen_patch_route_metadata(patch)
    denominator_data = hasproperty(patch, :denominator_data) ?
        copy(collect(patch.denominator_data)) :
        hasproperty(patch, :denominator_candidate) ?
        copy(collect(patch.denominator_candidate.raw_denominators)) :
        Any[]
    normalized_contribution_count = hasproperty(patch, :normalized_local_contributions) ?
        length(patch.normalized_local_contributions) :
        hasproperty(patch, :sequence_expansions) ?
        length(patch.sequence_expansions) :
        0
    return (;
        patch_size = patch.size,
        substitution_variable = patch.substitution_variable,
        denominator_data = denominator_data,
        local_certificate_count = hasproperty(patch, :local_certificates) ?
            length(patch.local_certificates) :
            length(patch.local_contributions),
        normalized_contribution_count = normalized_contribution_count,
        patch_replay_metadata = hasproperty(patch, :replay_metadata) ?
            patch.replay_metadata :
            nothing,
    )
end

function _polynomial_quillen_patch_route_core_verification(adapter)
    patch_verified_ok = verify_quillen_patch(adapter.quillen_patch)
    route_tag_ok = adapter.route == :quillen_patch

    recomputed_target_matrix = _polynomial_quillen_route_target_matrix(
        adapter.target,
        adapter.quillen_patch,
    )
    patch_target_matrix = _polynomial_quillen_patch_target_matrix(adapter.quillen_patch)
    recomputed_factors = copy(collect(_polynomial_quillen_patch_factors(adapter.quillen_patch)))
    recomputed_product = _polynomial_route_factor_product(
        recomputed_factors,
        base_ring(recomputed_target_matrix),
        nrows(recomputed_target_matrix),
    )
    replay_metadata = _polynomial_quillen_patch_route_metadata(adapter.quillen_patch)

    target_matrix_ok = recomputed_target_matrix == adapter.target_matrix
    patch_target_ok = patch_target_matrix == adapter.target_matrix
    factors_ok = _polynomial_route_factor_sequences_equal(
        adapter.global_elementary_factors,
        recomputed_factors,
    )
    product_ok =
        recomputed_product == adapter.product &&
        recomputed_product == adapter.target_matrix &&
        _polynomial_quillen_patch_product(adapter.quillen_patch) == recomputed_product
    replay_metadata_ok = replay_metadata == adapter.replay_metadata
    overall_core_ok =
        route_tag_ok &&
        patch_verified_ok &&
        target_matrix_ok &&
        patch_target_ok &&
        factors_ok &&
        product_ok &&
        replay_metadata_ok

    return (;
        route_tag_ok,
        patch_verified_ok,
        target_matrix_ok,
        patch_target_ok,
        factors_ok,
        product_ok,
        replay_metadata_ok,
        overall_core_ok,
    )
end

function _polynomial_quillen_patch_route_verification(adapter)
    core = _polynomial_quillen_patch_route_core_verification(adapter)
    stored_verification_ok = adapter.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function _verify_polynomial_quillen_patch_route_adapter(adapter)::Bool
    try
        verification = _polynomial_quillen_patch_route_verification(adapter)
        return verification.overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_quillen_patch_route_adapter(target, patch)
    verify_quillen_patch(patch) ||
        throw(ArgumentError("Quillen patch must verify before route adaptation"))
    target_matrix = _polynomial_quillen_route_target_matrix(target, patch)
    patch_target = _polynomial_quillen_patch_target_matrix(patch)
    patch_target == target_matrix ||
        throw(ArgumentError("Quillen patch target does not match route target"))
    factors = copy(collect(_polynomial_quillen_patch_factors(patch)))
    product = _polynomial_route_factor_product(
        factors,
        base_ring(target_matrix),
        nrows(target_matrix),
    )
    product == target_matrix ||
        throw(ArgumentError("Quillen patch route factors do not multiply to the target"))
    _polynomial_quillen_patch_product(patch) == product ||
        throw(ArgumentError("Quillen patch stored product does not match the adapted product"))
    replay_metadata = _polynomial_quillen_patch_route_metadata(patch)
    raw = PolynomialQuillenPatchRouteAdapter(
        target,
        :quillen_patch,
        patch,
        factors,
        product,
        target_matrix,
        replay_metadata,
        nothing,
    )
    verification = _polynomial_quillen_patch_route_core_verification(raw)
    adapter = PolynomialQuillenPatchRouteAdapter(
        target,
        :quillen_patch,
        patch,
        factors,
        product,
        target_matrix,
        replay_metadata,
        verification,
    )
    _verify_polynomial_quillen_patch_route_adapter(adapter) ||
        error("internal Quillen patch route adapter verification failed")
    return adapter
end

function _polynomial_staged_failure_evidence(A; allow_recursive_column_peel::Bool = true)
    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    if !_polynomial_exact_field_backed_ring(R)
        return (;
            error_type = :ArgumentError,
            message = _polynomial_unsupported_coefficient_ring_message(),
        )
    end

    X = _supported_local_sl3_generator(A, R, ring_profile)
    if X !== nothing
        return (; error_type = :none, message = "")
    end

    if _polynomial_quillen_supplied_evidence_route_certificate(A) !== nothing
        return (; error_type = :none, message = "")
    end

    if nrows(A) > 3
        staged_error = nothing
        try
            reduce_sln_to_sl3(A)
            return (; error_type = :none, message = "")
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            staged_error = err
        end

        if allow_recursive_column_peel
            try
                _polynomial_recursive_column_peel_route_certificate(
                    A;
                    route_tag = :polynomial_column_peel,
                )
                return (; error_type = :none, message = "")
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
            end
        end

        return (;
            error_type = Symbol(nameof(typeof(staged_error))),
            message = sprint(showerror, staged_error),
        )
    end

    try
        _throw_staged_factorization_failure(A, :polynomial, nothing)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return (;
            error_type = Symbol(nameof(typeof(err))),
            message = sprint(showerror, err),
        )
    end
end

function _polynomial_route_certificate(A, route::Symbol, factors, product, evidence, status::Symbol)
    stored_factors = copy(collect(factors))
    raw = PolynomialFactorizationRouteCertificate(
        A,
        route,
        stored_factors,
        product,
        evidence,
        status,
        nothing,
    )
    verification = _polynomial_factorization_route_core_verification(raw)
    certificate = PolynomialFactorizationRouteCertificate(
        A,
        route,
        stored_factors,
        product,
        evidence,
        status,
        verification,
    )
    _verify_polynomial_factorization_route_certificate(certificate) ||
        error("internal polynomial factorization route certificate verification failed")
    return certificate
end

function _verify_polynomial_factorization_route_certificate(cert)::Bool
    try
        verification = _polynomial_factorization_route_verification(cert)
        return verification.overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_factorization_route_verification(cert)
    core = _polynomial_factorization_route_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function _polynomial_factorization_route_core_verification(cert)
    A = cert.matrix
    route = cert.route
    route_tag_ok = route in _POLYNOMIAL_FACTORIZATION_ROUTE_TAGS
    square_ok = nrows(A) == ncols(A)
    n = square_ok ? nrows(A) : 0
    R = square_ok ? base_ring(A) : nothing

    ring_profile_ok = false
    determinant_ok = false
    if square_ok
        ring_profile_ok = _factorization_ring_profile(R) == :polynomial
        determinant_ok = ring_profile_ok && det(A) == one(R)
    end

    successful_route = route in (:fast_local_sl3, :disjoint_local_blocks, :quillen_patch) ||
        _is_polynomial_column_peel_route(route)
    supported_status_ok = successful_route && cert.status == :supported
    staged_status_ok = route == :staged_failure && cert.status == :staged
    status_ok = supported_status_ok || staged_status_ok
    factors_vector_ok = cert.factors isa AbstractVector

    replayed_product =
        square_ok && ring_profile_ok && factors_vector_ok ?
        _polynomial_route_factor_product(cert.factors, R, n) :
        nothing
    product_replay_ok = replayed_product !== nothing
    product_matches_stored_ok = product_replay_ok && replayed_product == cert.product
    product_matches_matrix_ok = successful_route && product_replay_ok && replayed_product == A
    factorization_ok = successful_route && factors_vector_ok && verify_factorization(A, cert.factors)
    staged_product_ok =
        route == :staged_failure &&
        square_ok &&
        cert.product == identity_matrix(R, n) &&
        isempty(cert.factors)
    evidence_ok = _polynomial_route_evidence_ok(cert)

    successful_route_ok =
        successful_route &&
        supported_status_ok &&
        product_matches_stored_ok &&
        product_matches_matrix_ok &&
        factorization_ok &&
        evidence_ok
    staged_route_ok =
        route == :staged_failure &&
        staged_status_ok &&
        staged_product_ok &&
        evidence_ok
    overall_core_ok =
        route_tag_ok &&
        square_ok &&
        ring_profile_ok &&
        determinant_ok &&
        status_ok &&
        (successful_route_ok || staged_route_ok)

    return (;
        route_tag_ok,
        square_ok,
        ring_profile_ok,
        determinant_ok,
        status_ok,
        factors_vector_ok,
        product_replay_ok,
        product_matches_stored_ok,
        product_matches_matrix_ok,
        factorization_ok,
        staged_product_ok,
        evidence_ok,
        successful_route_ok,
        staged_route_ok,
        overall_core_ok,
    )
end

function _polynomial_route_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        nrows(factor) == n || throw(ArgumentError("route certificate factor has wrong row count"))
        ncols(factor) == n || throw(ArgumentError("route certificate factor has wrong column count"))
        _same_base_ring(base_ring(factor), R) ||
            throw(ArgumentError("route certificate factor has wrong base ring"))
        product *= factor
    end
    return product
end

function _polynomial_route_factor_sequences_equal(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left, right)
        left[idx] == right[idx] || return false
    end
    return true
end

function _polynomial_route_evidence_ok(cert)::Bool
    try
        if cert.route == :fast_local_sl3
            return cert.evidence isa SL3LocalRealizationCertificate &&
                cert.evidence.target == cert.matrix &&
                verify_sl3_local_realization(cert.evidence) &&
                _polynomial_route_factor_sequences_equal(cert.factors, cert.evidence.factors)
        elseif cert.route == :disjoint_local_blocks
            return cert.evidence isa SLNToSL3Reduction &&
                cert.evidence.original_matrix == cert.matrix &&
                cert.evidence.product == cert.matrix &&
                verify_sln_to_sl3_reduction(cert.evidence) &&
                _polynomial_route_factor_sequences_equal(cert.factors, cert.evidence.factors)
        elseif _is_polynomial_column_peel_route(cert.route)
            return cert.evidence isa PolynomialColumnPeelCertificate &&
                cert.evidence.original_matrix == cert.matrix &&
                cert.evidence.product == cert.matrix &&
                _verify_polynomial_column_peel_certificate(cert.evidence) &&
                _polynomial_route_factor_sequences_equal(cert.factors, cert.evidence.factors)
        elseif cert.route == :quillen_patch
            return cert.evidence isa PolynomialQuillenPatchRouteAdapter &&
                cert.evidence.target_matrix == cert.matrix &&
                _verify_polynomial_quillen_patch_route_adapter(cert.evidence) &&
                _polynomial_route_factor_sequences_equal(
                    cert.factors,
                    cert.evidence.global_elementary_factors,
                )
        elseif cert.route == :staged_failure
            hasproperty(cert.evidence, :error_type) &&
                hasproperty(cert.evidence, :message) &&
                cert.evidence.error_type isa Symbol &&
                cert.evidence.message isa AbstractString &&
                !isempty(cert.evidence.message) ||
                return false
            fresh_evidence = _polynomial_staged_failure_evidence(cert.matrix)
            return cert.evidence == fresh_evidence &&
                fresh_evidence.error_type == :ArgumentError &&
                !isempty(fresh_evidence.message)
        end
    catch err
        err isa InterruptException && rethrow()
        return false
    end

    return false
end

function _polynomial_verified_route_factors(A)
    certificate = _polynomial_factorization_route_certificate(A)
    return _polynomial_verified_certificate_factors(certificate)
end

function _polynomial_verified_certificate_factors(certificate)
    if certificate.status == :supported
        if _verify_polynomial_factorization_route_certificate(certificate) &&
                verify_factorization(certificate.matrix, certificate.factors)
            return certificate.factors
        end
        error("internal polynomial factorization route certificate verification failed")
    elseif certificate.status == :staged
        _throw_polynomial_staged_certificate_failure(certificate)
    end

    throw(ArgumentError("unsupported polynomial factorization route certificate status $(certificate.status)"))
end

function _throw_polynomial_staged_certificate_failure(certificate)
    evidence = certificate.evidence
    if hasproperty(evidence, :message) &&
            evidence.message isa AbstractString &&
            !isempty(evidence.message)
        if occursin("ordinary polynomial reduction currently requires a univariate base ring", evidence.message)
            _throw_staged_factorization_failure(certificate.matrix, :polynomial, nothing)
        end
        throw(ArgumentError(evidence.message))
    end

    _throw_staged_factorization_failure(certificate.matrix, :polynomial, nothing)
end

function _laurent_sl_fallback_factorization(A)
    try
        reduction = reduce_sln_to_sl3(A)
        verify_factorization(A, reduction.factors) && return reduction.factors
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
    end

    certificate = _factor_laurent_sl_column_peel(A)
    verify_factorization(A, certificate.factors) && return certificate.factors
    error("internal Laurent column-peel factorization failed exact verification")
end

function elementary_factorization(A)
    _validate_factorization_matrix(A)

    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    normalized_A, normalization = _normalize_factorization_input(A, ring_profile)

    if ring_profile == :polynomial
        _require_polynomial_sl_determinant(normalized_A)
        return _polynomial_verified_route_factors(normalized_A)
    end

    if normalization.determinant_classification == :one
        return _laurent_sl_fallback_factorization(A)
    end

    _throw_staged_factorization_failure(normalized_A, ring_profile, normalization)
end

function verify_factorization(A, factors)::Bool
    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))

    R = base_ring(A)
    product = identity_matrix(R, nrows(A))

    for factor in factors
        nrows(factor) == nrows(A) || throw(ArgumentError("all factors must have the same size as A"))
        ncols(factor) == ncols(A) || throw(ArgumentError("all factors must have the same size as A"))
        _same_base_ring(base_ring(factor), R) || throw(ArgumentError("all factors must lie in the same base ring as A"))
        product *= factor
    end

    return product == A
end
