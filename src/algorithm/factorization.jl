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
    creation_identity::NamedTuple
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

struct SL3LocalFormWitnessSelection
    context::SL3RealizationInputContext
    selected_variable
    selected_variable_index
    selected_variable_name
    entries
    local_form_matrix
    monicity_witness::NamedTuple
    local_form_witness
    variable_change_metadata
    variable_change_status::Symbol
    normality_conjugation_metadata
    normality_conjugation_status::Symbol
    replay_status::Symbol
    support_status::Symbol
    witness_source::Symbol
    staged_diagnostic::NamedTuple
    verification
end

struct SL3MurthyQuillenLocalEvidenceProvider
    context::SL3RealizationInputContext
    witness_selection::SL3LocalFormWitnessSelection
    murthy_context
    murthy_certificate
    murthy_adapter
    quillen_local_sequences::Vector
    selected_variable
    selected_variable_index
    local_product
    denominator_metadata
    witness_metadata
    replay_metadata
    staged_diagnostic::NamedTuple
    metadata
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

struct PolynomialSL3QuillenMurthyRouteEvidence
    target
    route::Symbol
    context::SL3RealizationInputContext
    witness_selection::SL3LocalFormWitnessSelection
    local_evidence_provider::SL3MurthyQuillenLocalEvidenceProvider
    quillen_consumption
    quillen_route_adapter::PolynomialQuillenPatchRouteAdapter
    base_term_policy::Symbol
    base_term_factors::Vector
    replay_metadata
    verification
end

struct PolynomialSL3SuppliedQuillenRouteEvidence
    target
    route::Symbol
    supplied_evidence
    quillen_route_adapter::PolynomialQuillenPatchRouteAdapter
    replay_metadata::NamedTuple
    verification
end

struct PolynomialSL3IdentityQuillenRouteEvidence
    target
    route::Symbol
    factors::Vector
    product
    replay_metadata::NamedTuple
    verification
end

struct PolynomialSL3QuillenMurthyRouteConsumption
    original_input
    ring
    size::Int
    selected_variable
    murthy_adapters::Vector
    local_sequence_certificates::Vector
    patch
    raw_consumption
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

function _sl3_realization_input_context_matrix_snapshot(matrix)
    return (;
        ring = base_ring(matrix),
        rows = nrows(matrix),
        cols = ncols(matrix),
        entries = Tuple(matrix[row, col] for row in 1:nrows(matrix), col in 1:ncols(matrix)),
    )
end

function _sl3_realization_input_context_creation_snapshot(
    matrix,
    selected_variable,
    catalog_metadata,
    local_form_witness,
)
    return (;
        matrix = _sl3_realization_input_context_matrix_snapshot(matrix),
        selected_variable,
        catalog_metadata = deepcopy(catalog_metadata),
        local_form_witness = deepcopy(local_form_witness),
    )
end

function _sl3_realization_input_context_creation_identity(
    matrix,
    selected_variable,
    catalog_metadata,
    local_form_witness,
)
    return (;
        source = :sl3_realization_input_context,
        snapshot = _sl3_realization_input_context_creation_snapshot(
            matrix,
            selected_variable,
            catalog_metadata,
            local_form_witness,
        ),
    )
end

function _sl3_realization_input_context_creation_identity_ok(
    creation_identity,
    matrix,
    selected_variable,
    catalog_metadata,
    local_form_witness,
)::Bool
    creation_identity isa NamedTuple || return false
    hasproperty(creation_identity, :source) || return false
    creation_identity.source == :sl3_realization_input_context || return false
    hasproperty(creation_identity, :snapshot) || return false
    return creation_identity.snapshot == _sl3_realization_input_context_creation_snapshot(
        matrix,
        selected_variable,
        catalog_metadata,
        local_form_witness,
    )
end

function _sl3_realization_input_context_fields(
    A;
    selected_variable = nothing,
    catalog_metadata = (;),
    creation_identity = nothing,
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
        creation_identity = creation_identity === nothing ?
            _sl3_realization_input_context_creation_identity(
                A,
                selected,
                catalog_metadata,
                local_form_witness,
            ) :
            creation_identity,
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
        creation_identity = context.creation_identity,
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
    creation_identity_ok =
        context.creation_identity == recomputed.creation_identity &&
        _sl3_realization_input_context_creation_identity_ok(
            context.creation_identity,
            context.matrix,
            context.selected_variable,
            context.catalog_metadata,
            context.local_form_witness,
        )
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
        creation_identity_ok &&
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
        creation_identity_ok,
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

function _sl3_local_witness_hint_or_context(hint, stored)
    return hint === nothing ? stored : hint
end

function _sl3_local_witness_selected_variable(context, hint)
    R = context.base_ring
    selected_hint = _sl3_local_witness_hint_or_context(hint, context.selected_variable)
    selected, selected_index, status =
        _sl3_realization_input_context_selected_variable(R, selected_hint)
    status == :passes ||
        throw(ArgumentError("SL_3 local witness selection requires a selected variable"))
    if context.selected_variable !== nothing && selected != context.selected_variable
        throw(ArgumentError("selected variable hint does not match the SL_3 context"))
    end
    if context.selected_variable_index !== nothing &&
            selected_index != context.selected_variable_index
        throw(ArgumentError("selected variable index does not match the SL_3 context"))
    end
    selected_name = string(collect(gens(R))[selected_index])
    return selected, selected_index, selected_name
end

function _sl3_local_witness_entries_from_data(data)
    data === nothing && return nothing
    entries = _sl3_realization_input_context_extract(data, (:entries, :local_form_entries))
    entries !== nothing && return entries
    if all(field -> hasproperty(data, field), (:p, :q, :r, :s))
        return (; p = data.p, q = data.q, r = data.r, s = data.s)
    end
    return nothing
end

function _sl3_local_witness_matrix_from_entries(R, entries)
    entries === nothing && return nothing
    return matrix(R, [
        entries.p entries.q zero(R);
        entries.r entries.s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_local_witness_extract_matrix(data)
    return _sl3_realization_input_context_extract(
        data,
        (:local_form_matrix, :special_form_matrix, :transformed_matrix, :local_target),
    )
end

function _sl3_local_witness_source_matrix(data)
    return _sl3_realization_input_context_extract(
        data,
        (:source_matrix, :input_matrix, :context_matrix, :original_matrix),
    )
end

function _sl3_local_witness_same_ring_entries(R, entries, source_label)
    for key in (:p, :q, :r, :s)
        hasproperty(entries, key) ||
            throw(ArgumentError("$(source_label) local-form entries are malformed"))
        parent(getproperty(entries, key)) == R ||
            throw(ArgumentError("$(source_label) local-form entries must lie in the SL_3 context ring"))
    end
    return entries
end

function _sl3_local_witness_checked_entries(
    R,
    matrix_value,
    selected_variable,
    selected_index,
    source_label,
)
    base_ring(matrix_value) == R ||
        throw(ArgumentError("$(source_label) local-form target must lie in the SL_3 context ring"))
    nrows(matrix_value) == 3 && ncols(matrix_value) == 3 ||
        throw(ArgumentError("$(source_label) local-form target must be a 3 x 3 special-form matrix"))
    entries = _sl3_local_target_entries(matrix_value)
    entries !== nothing ||
        throw(ArgumentError("$(source_label) local-form target is not in SL_3 special form"))
    det(matrix_value) == one(R) ||
        throw(ArgumentError("$(source_label) local-form target must have determinant 1"))
    monicity_witness = _sl3_local_monicity_witness(entries.p, selected_index, R)
    monicity_witness.is_monic ||
        throw(ArgumentError("$(source_label) local-form p is not monic in the selected variable"))
    monicity_witness.variable == selected_variable ||
        throw(ArgumentError("$(source_label) monicity witness variable does not match the selected variable"))
    return entries, matrix_value, monicity_witness
end

function _sl3_local_witness_status_from_data(
    R,
    data,
    selected_variable,
    selected_index,
    source_label,
)
    entries = _sl3_local_witness_entries_from_data(data)
    matrix_value = _sl3_local_witness_extract_matrix(data)

    entries === nothing && matrix_value === nothing &&
        return nothing, nothing, nothing

    checked_entries = nothing
    checked_matrix = nothing
    checked_monicity = nothing

    if matrix_value !== nothing
        checked_entries, checked_matrix, checked_monicity =
            _sl3_local_witness_checked_entries(
                R,
                matrix_value,
                selected_variable,
                selected_index,
                source_label,
            )
    end

    if entries !== nothing
        checked_entry_data = _sl3_local_witness_same_ring_entries(R, entries, source_label)
        entry_matrix = _sl3_local_witness_matrix_from_entries(R, checked_entry_data)
        entry_entries, entry_matrix, entry_monicity =
            _sl3_local_witness_checked_entries(
                R,
                entry_matrix,
                selected_variable,
                selected_index,
                source_label,
            )
        if checked_matrix !== nothing
            checked_matrix == entry_matrix && checked_entries == entry_entries ||
                throw(ArgumentError("$(source_label) local-form entries do not match the supplied local-form matrix"))
        else
            checked_entries = entry_entries
            checked_matrix = entry_matrix
            checked_monicity = entry_monicity
        end
    end

    return checked_entries, checked_matrix, checked_monicity
end

function _sl3_local_witness_metadata_status(
    context,
    metadata,
    selected_variable,
    selected_index,
    source_label,
)
    metadata === nothing && return :missing, nothing, nothing, nothing

    source_matrix = _sl3_local_witness_source_matrix(metadata)
    if source_matrix !== nothing && source_matrix != context.matrix
        throw(ArgumentError("$(source_label) source matrix does not match the SL_3 context"))
    end

    selected_variable_hint = _sl3_realization_input_context_extract(
        metadata,
        (:selected_variable, :selected_generator, :generator, :variable),
    )
    if selected_variable_hint !== nothing
        replay_selected, replay_index, replay_status =
            _sl3_realization_input_context_selected_variable(context.base_ring, selected_variable_hint)
        replay_status == :passes ||
            throw(ArgumentError("$(source_label) selected variable is not valid for the SL_3 context"))
        replay_selected == selected_variable && replay_index == selected_index ||
            throw(ArgumentError("$(source_label) selected variable does not match the SL_3 context"))
    end

    entries, matrix_value, monicity_witness = _sl3_local_witness_status_from_data(
        context.base_ring,
        metadata,
        selected_variable,
        selected_index,
        source_label,
    )

    has_replay_payload = _sl3_realization_input_context_has_replay_payload(metadata)
    has_context_matrix = source_matrix !== nothing && source_matrix == context.matrix
    if entries === nothing
        return :recorded, nothing, nothing, nothing
    end

    if has_replay_payload && has_context_matrix
        selected_variable_hint !== nothing ||
            return :recorded, entries, matrix_value, monicity_witness
        return :replayed, entries, matrix_value, monicity_witness
    end

    return :recorded, entries, matrix_value, monicity_witness
end

function _sl3_local_witness_replay_matches(
    entries,
    matrix_value,
    monicity_witness,
    reference_entries,
    reference_matrix,
    reference_monicity,
)
    return entries == reference_entries &&
        matrix_value == reference_matrix &&
        monicity_witness == reference_monicity
end

function _sl3_local_witness_require_replay_consistency(
    status::Symbol,
    entries,
    matrix_value,
    monicity_witness,
    reference_entries,
    reference_matrix,
    reference_monicity,
    source_label,
)
    status == :replayed || return nothing
    _sl3_local_witness_replay_matches(
        entries,
        matrix_value,
        monicity_witness,
        reference_entries,
        reference_matrix,
        reference_monicity,
    ) || throw(ArgumentError("$(source_label) replayed local-form witness conflicts with another replayed witness"))
    return nothing
end

function _sl3_local_form_witness_selection_staged_diagnostic(
    variable_change_status::Symbol,
    normality_conjugation_status::Symbol,
    support_status::Symbol,
)
    missing_evidence = [:local_form]
    partial_evidence = Symbol[]

    for (label, status) in (
        (:variable_change, variable_change_status),
        (:normality_conjugation, normality_conjugation_status),
    )
        if status == :recorded
            push!(partial_evidence, label)
        elseif status == :missing
            push!(missing_evidence, label)
        end
    end

    reason = support_status == :supported ?
        "supported local-form witness selected" :
        "missing supported local-form witness"
    message = support_status == :supported ?
        "SL_3 local-form witness is replayable" :
        "SL_3 local-form witness is staged"

    return (;
        reason,
        message,
        missing_evidence = Tuple(missing_evidence),
        partial_evidence = Tuple(partial_evidence),
        status = support_status,
    )
end

function _sl3_local_form_witness_selection_fields(
    context::SL3RealizationInputContext;
    selected_variable = nothing,
    local_form_witness = nothing,
    variable_change_metadata = nothing,
    normality_conjugation_metadata = nothing,
)
    _verify_sl3_realization_input_context(context) ||
        throw(ArgumentError("SL_3 local witness selection requires a verified realization input context"))

    selected, selected_index, selected_name =
        _sl3_local_witness_selected_variable(context, selected_variable)
    local_form_data = _sl3_local_witness_hint_or_context(local_form_witness, context.local_form_witness)
    variable_change_data =
        _sl3_local_witness_hint_or_context(variable_change_metadata, context.variable_change_metadata)
    normality_data = _sl3_local_witness_hint_or_context(
        normality_conjugation_metadata,
        context.normality_conjugation_metadata,
    )

    context_entries = _sl3_local_target_entries(context.matrix)
    context_supported = false
    context_monicity = nothing
    if context_entries !== nothing
        context_entries, context_matrix, context_monicity = _sl3_local_witness_checked_entries(
            context.base_ring,
            context.matrix,
            selected,
            selected_index,
            "context matrix",
        )
        context_supported = true

        witness_entries, witness_matrix, _ = _sl3_local_witness_status_from_data(
            context.base_ring,
            local_form_data,
            selected,
            selected_index,
            "supplied local-form witness",
        )
        if witness_entries !== nothing
            witness_entries == context_entries && witness_matrix == context_matrix ||
                throw(ArgumentError("supplied local-form witness does not match the context special form"))
        end
    else
        _sl3_local_witness_status_from_data(
            context.base_ring,
            local_form_data,
            selected,
            selected_index,
            "supplied local-form witness",
        )
    end

    variable_change_status, variable_entries, variable_matrix, variable_monicity =
        _sl3_local_witness_metadata_status(
            context,
            variable_change_data,
            selected,
            selected_index,
            "variable-change metadata",
        )
    normality_status, normality_entries, normality_matrix, normality_monicity =
        _sl3_local_witness_metadata_status(
            context,
            normality_data,
            selected,
            selected_index,
            "normality/conjugation metadata",
        )

    if context_supported
        _sl3_local_witness_require_replay_consistency(
            variable_change_status,
            variable_entries,
            variable_matrix,
            variable_monicity,
            context_entries,
            context_matrix,
            context_monicity,
            "variable-change metadata",
        )
        _sl3_local_witness_require_replay_consistency(
            normality_status,
            normality_entries,
            normality_matrix,
            normality_monicity,
            context_entries,
            context_matrix,
            context_monicity,
            "normality/conjugation metadata",
        )
    elseif variable_change_status == :replayed && normality_status == :replayed
        _sl3_local_witness_require_replay_consistency(
            normality_status,
            normality_entries,
            normality_matrix,
            normality_monicity,
            variable_entries,
            variable_matrix,
            variable_monicity,
            "normality/conjugation metadata",
        )
    end

    replay_status = :missing
    support_status = :staged
    witness_source = :staged
    entries = nothing
    local_form_matrix = nothing
    monicity_witness = (
        variable = selected,
        variable_index = selected_index,
        degree = -1,
        leading_coefficient = zero(context.base_ring),
        is_monic = false,
    )

    if context_supported
        replay_status = :replayed
        support_status = :supported
        witness_source = :already_special_form
        entries = context_entries
        local_form_matrix = context.matrix
        monicity_witness = context_monicity
    elseif variable_change_status == :replayed
        replay_status = :replayed
        support_status = :supported
        witness_source = :variable_change
        entries = variable_entries
        local_form_matrix = variable_matrix
        monicity_witness = variable_monicity
    elseif normality_status == :replayed
        replay_status = :replayed
        support_status = :supported
        witness_source = :normality_conjugation
        entries = normality_entries
        local_form_matrix = normality_matrix
        monicity_witness = normality_monicity
    end

    staged_diagnostic = _sl3_local_form_witness_selection_staged_diagnostic(
        variable_change_status,
        normality_status,
        support_status,
    )

    return (;
        context,
        selected_variable = selected,
        selected_variable_index = selected_index,
        selected_variable_name = selected_name,
        entries,
        local_form_matrix,
        monicity_witness,
        local_form_witness = local_form_data,
        variable_change_metadata = variable_change_data,
        variable_change_status,
        normality_conjugation_metadata = normality_data,
        normality_conjugation_status = normality_status,
        replay_status,
        support_status,
        witness_source,
        staged_diagnostic,
    )
end

function _sl3_local_form_witness_selection_core_verification(selection)
    recomputed = _sl3_local_form_witness_selection_fields(
        selection.context;
        selected_variable = selection.selected_variable,
        local_form_witness = selection.local_form_witness,
        variable_change_metadata = selection.variable_change_metadata,
        normality_conjugation_metadata = selection.normality_conjugation_metadata,
    )

    context_ok = selection.context == recomputed.context
    selected_variable_ok = selection.selected_variable == recomputed.selected_variable
    selected_variable_index_ok =
        selection.selected_variable_index == recomputed.selected_variable_index
    selected_variable_name_ok =
        selection.selected_variable_name == recomputed.selected_variable_name
    entries_ok = selection.entries == recomputed.entries
    local_form_matrix_ok = selection.local_form_matrix == recomputed.local_form_matrix
    monicity_witness_ok = selection.monicity_witness == recomputed.monicity_witness
    local_form_witness_ok = selection.local_form_witness == recomputed.local_form_witness
    variable_change_metadata_ok =
        selection.variable_change_metadata == recomputed.variable_change_metadata
    variable_change_status_ok =
        selection.variable_change_status == recomputed.variable_change_status
    normality_conjugation_metadata_ok =
        selection.normality_conjugation_metadata == recomputed.normality_conjugation_metadata
    normality_conjugation_status_ok =
        selection.normality_conjugation_status == recomputed.normality_conjugation_status
    replay_status_ok = selection.replay_status == recomputed.replay_status
    support_status_ok = selection.support_status == recomputed.support_status
    witness_source_ok = selection.witness_source == recomputed.witness_source
    staged_diagnostic_ok = selection.staged_diagnostic == recomputed.staged_diagnostic

    overall_core_ok =
        context_ok &&
        selected_variable_ok &&
        selected_variable_index_ok &&
        selected_variable_name_ok &&
        entries_ok &&
        local_form_matrix_ok &&
        monicity_witness_ok &&
        local_form_witness_ok &&
        variable_change_metadata_ok &&
        variable_change_status_ok &&
        normality_conjugation_metadata_ok &&
        normality_conjugation_status_ok &&
        replay_status_ok &&
        support_status_ok &&
        witness_source_ok &&
        staged_diagnostic_ok

    return (;
        context_ok,
        selected_variable_ok,
        selected_variable_index_ok,
        selected_variable_name_ok,
        entries_ok,
        local_form_matrix_ok,
        monicity_witness_ok,
        local_form_witness_ok,
        variable_change_metadata_ok,
        variable_change_status_ok,
        normality_conjugation_metadata_ok,
        normality_conjugation_status_ok,
        replay_status_ok,
        support_status_ok,
        witness_source_ok,
        staged_diagnostic_ok,
        overall_core_ok,
    )
end

function _sl3_local_form_witness_selection_verification(selection)
    core = _sl3_local_form_witness_selection_core_verification(selection)
    stored_verification_ok = selection.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _select_sl3_local_form_witness(
    context::SL3RealizationInputContext;
    selected_variable = nothing,
    local_form_witness = nothing,
    variable_change_metadata = nothing,
    normality_conjugation_metadata = nothing,
)
    fields = _sl3_local_form_witness_selection_fields(
        context;
        selected_variable,
        local_form_witness,
        variable_change_metadata,
        normality_conjugation_metadata,
    )
    unchecked = SL3LocalFormWitnessSelection(values(merge(fields, (; verification = nothing,)))...)
    verification = _sl3_local_form_witness_selection_core_verification(unchecked)
    checked = SL3LocalFormWitnessSelection(values(merge(fields, (; verification,)))...)
    _verify_sl3_local_form_witness_selection(checked) ||
        error("internal SL_3 local witness selection verification failed")
    return checked
end

function _verify_sl3_local_form_witness_selection(selection)::Bool
    try
        return _sl3_local_form_witness_selection_verification(selection).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sl3_murthy_quillen_local_witness_metadata(
    selection::SL3LocalFormWitnessSelection,
)
    return (;
        context_metadata = selection.context.catalog_metadata,
        local_form_witness = selection.local_form_witness,
        variable_change_metadata = selection.variable_change_metadata,
        normality_conjugation_metadata = selection.normality_conjugation_metadata,
        witness_source = selection.witness_source,
        replay_status = selection.replay_status,
        support_status = selection.support_status,
        selected_variable = selection.selected_variable,
        selected_variable_index = selection.selected_variable_index,
        selected_variable_name = selection.selected_variable_name,
    )
end

function _sl3_murthy_quillen_local_provider_inputs(
    ;
    witness,
    local_unit_witnesses,
    split_witness,
    bezout_witness,
)
    return (;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
end

function _sl3_murthy_quillen_local_witness_payload_ok(
    selection::SL3LocalFormWitnessSelection,
)::Bool
    data = selection.local_form_witness
    data === nothing && return false
    _sl3_realization_input_context_has_replay_payload(data) || return false

    source_matrix = _sl3_local_witness_source_matrix(data)
    source_matrix !== nothing && source_matrix == selection.context.matrix || return false

    selected_variable_hint = _sl3_realization_input_context_extract(
        data,
        (:selected_variable, :selected_generator, :generator, :variable),
    )
    selected_variable_hint !== nothing || return false
    replay_selected, replay_index, replay_status =
        _sl3_realization_input_context_selected_variable(
            selection.context.base_ring,
            selected_variable_hint,
        )
    replay_status == :passes || return false
    replay_selected == selection.selected_variable || return false
    replay_index == selection.selected_variable_index || return false

    entries, local_form_matrix, monicity_witness = _sl3_local_witness_status_from_data(
        selection.context.base_ring,
        data,
        selection.selected_variable,
        selection.selected_variable_index,
        "SL_3 Murthy/Quillen local evidence provider #236 witness",
    )
    return entries == selection.entries &&
           local_form_matrix == selection.local_form_matrix &&
           monicity_witness == selection.monicity_witness
end

function _sl3_murthy_quillen_local_denominator_metadata(
    adapter,
    quillen_local_sequences,
)
    factor_denominators = [
        factor.denominator for factor in adapter.local_factor_replay.factors
    ]
    sequence_denominator_products = if adapter.mode == :ordinary_quillen_factor_sequence
        [sequence.product_denominator for sequence in quillen_local_sequences]
    else
        Any[]
    end
    return (;
        factor_denominators,
        denominator_product = adapter.local_factor_replay.denominator_product,
        adapter_mode = adapter.mode,
        sequence_denominator_products,
    )
end

function _sl3_murthy_quillen_local_staged_diagnostic(
    adapter,
    quillen_local_sequences,
)
    if adapter.mode == :ordinary_quillen_factor_sequence
        return (;
            reason = :ordinary_quillen_sequences,
            message = "Murthy local SL_3 evidence provider produced verified ordinary Quillen local sequences",
            adapter_mode = adapter.mode,
            quillen_local_sequence_count = length(quillen_local_sequences),
            status = :supported,
        )
    end

    return (;
        reason = :localized_denominator_cleared_handoff,
        message = "Murthy local SL_3 evidence provider captured localized denominator-cleared replay; ordinary Quillen local sequence materialization remains staged",
        adapter_mode = adapter.mode,
        quillen_local_sequence_count = length(quillen_local_sequences),
        status = :staged,
    )
end

function _sl3_murthy_quillen_local_replay_metadata(
    selection::SL3LocalFormWitnessSelection,
    provider_inputs,
    denominator_metadata,
    staged_diagnostic,
    metadata,
)
    return (;
        original_matrix = selection.context.matrix,
        local_product = selection.local_form_matrix,
        context_metadata = selection.context.catalog_metadata,
        witness_source = selection.witness_source,
        provider_inputs,
        denominator_metadata,
        staged_diagnostic,
        provider_metadata = metadata,
    )
end

function _sl3_murthy_quillen_local_evidence_provider_fields(
    selection::SL3LocalFormWitnessSelection;
    witness = nothing,
    local_unit_witnesses = (;),
    split_witness = nothing,
    bezout_witness = nothing,
    metadata = (;),
)
    _verify_sl3_local_form_witness_selection(selection) ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider requires a verified local-form witness selection"))
    selection.support_status == :supported ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider requires supported local-form evidence"))
    selection.replay_status == :replayed ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider requires replayable local-form evidence"))
    selection.entries !== nothing ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider requires selected special-form entries"))
    selection.local_form_matrix !== nothing ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider requires a selected local-form matrix"))
    selection.local_form_witness == selection.context.local_form_witness ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider requires #236 witness metadata bound to the realization context"))
    _sl3_murthy_quillen_local_witness_payload_ok(selection) ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider requires a replayable #236 local-form witness payload"))
    selection.local_form_matrix == selection.context.matrix ||
        throw(ArgumentError("SL_3 Murthy/Quillen local evidence provider currently requires the selected local product to match the original SL_3 context matrix; transformed local-form replay remains staged"))

    provider_inputs = _sl3_murthy_quillen_local_provider_inputs(
        ;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
    witness_metadata = _sl3_murthy_quillen_local_witness_metadata(selection)
    murthy_context = sl3_local_murthy_input_context(
        selection.entries.p,
        selection.entries.q,
        selection.entries.r,
        selection.entries.s,
        selection.selected_variable;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
    murthy_certificate = realize_sl3_local_certificate(murthy_context)
    murthy_adapter = _murthy_quillen_local_adapter(
        murthy_certificate,
        selection.local_form_matrix,
        selection.selected_variable;
        witness_metadata,
    )
    quillen_local_sequences = murthy_adapter.mode == :ordinary_quillen_factor_sequence ?
        quillen_local_sequences_from_murthy_adapters(
            selection.local_form_matrix,
            selection.selected_variable,
            [murthy_adapter],
        ) :
        Any[]
    denominator_metadata = _sl3_murthy_quillen_local_denominator_metadata(
        murthy_adapter,
        quillen_local_sequences,
    )
    staged_diagnostic = _sl3_murthy_quillen_local_staged_diagnostic(
        murthy_adapter,
        quillen_local_sequences,
    )
    replay_metadata = _sl3_murthy_quillen_local_replay_metadata(
        selection,
        provider_inputs,
        denominator_metadata,
        staged_diagnostic,
        metadata,
    )

    return (;
        context = selection.context,
        witness_selection = selection,
        murthy_context,
        murthy_certificate,
        murthy_adapter,
        quillen_local_sequences,
        selected_variable = selection.selected_variable,
        selected_variable_index = selection.selected_variable_index,
        local_product = selection.local_form_matrix,
        denominator_metadata,
        witness_metadata,
        replay_metadata,
        staged_diagnostic,
        metadata,
    )
end

function _same_cached_provenance_data(left, right)::Bool
    left === right && return true
    typeof(left) === typeof(right) || return false

    values_equal = try
        left == right
    catch err
        err isa InterruptException && rethrow()
        false
    end
    values_equal && return true

    if left isa NamedTuple
        keys(left) == keys(right) || return false
        for key in keys(left)
            _same_cached_provenance_data(
                getproperty(left, key),
                getproperty(right, key),
            ) || return false
        end
        return true
    elseif left isa Tuple
        length(left) == length(right) || return false
        for idx in eachindex(left, right)
            _same_cached_provenance_data(left[idx], right[idx]) || return false
        end
        return true
    elseif left isa AbstractArray
        size(left) == size(right) || return false
        for idx in eachindex(left, right)
            _same_cached_provenance_data(left[idx], right[idx]) || return false
        end
        return true
    elseif Base.isstructtype(typeof(left)) && fieldcount(typeof(left)) > 0
        for name in fieldnames(typeof(left))
            _same_cached_provenance_data(
                getfield(left, name),
                getfield(right, name),
            ) || return false
        end
        return true
    end

    return false
end

function _same_sl3_local_realization_certificate_data(left, right)::Bool
    return left.target == right.target &&
           left.branch == right.branch &&
           _polynomial_route_factor_sequences_equal(left.factors, right.factors) &&
           left.selected_variable == right.selected_variable &&
           _same_cached_provenance_data(left.witness, right.witness)
end

function _same_sl3_local_murthy_input_context_data(left, right)::Bool
    return left.R === right.R &&
           left.X == right.X &&
           left.var_idx == right.var_idx &&
           left.entries == right.entries &&
           left.target == right.target &&
           left.determinant == right.determinant &&
           left.degree_p == right.degree_p &&
           left.degree_q == right.degree_q &&
           left.p0 == right.p0 &&
           left.q0 == right.q0 &&
           left.p_monic == right.p_monic &&
           left.global_units == right.global_units &&
           left.local_units == right.local_units &&
           left.local_unit_witnesses == right.local_unit_witnesses &&
           left.split_witness == right.split_witness &&
           left.bezout_witness == right.bezout_witness
end

function _same_murthy_quillen_local_adapter_data(left, right)::Bool
    quillen_sequence_presence_ok =
        (left.quillen_factor_sequence === nothing) ==
        (right.quillen_factor_sequence === nothing)
    quillen_local_presence_ok =
        (left.quillen_local_certificate === nothing) ==
        (right.quillen_local_certificate === nothing)
    return left.original_input == right.original_input &&
           left.ring == right.ring &&
           left.size == right.size &&
           left.selected_variable == right.selected_variable &&
           _same_sl3_local_realization_certificate_data(
               left.murthy_certificate,
               right.murthy_certificate,
           ) &&
           _same_sl3_local_factor_replay(left.local_factor_replay, right.local_factor_replay) &&
           left.mode == right.mode &&
           left.materialized_factors == right.materialized_factors &&
           left.local_product == right.local_product &&
           left.local_correction == right.local_correction &&
           quillen_sequence_presence_ok &&
           quillen_local_presence_ok &&
           left.witness_metadata == right.witness_metadata &&
           left.replay_metadata == right.replay_metadata
end

function _sl3_murthy_quillen_local_evidence_provider_core_verification(provider)
    selection_ok = _verify_sl3_local_form_witness_selection(provider.witness_selection)
    context_ok =
        provider.context == provider.witness_selection.context &&
        _verify_sl3_realization_input_context(provider.context)

    provider_inputs =
        hasproperty(provider.replay_metadata, :provider_inputs) ?
        provider.replay_metadata.provider_inputs :
        hasproperty(provider.metadata, :provider_inputs) ?
        provider.metadata.provider_inputs :
        _sl3_murthy_quillen_local_provider_inputs(
            ;
            witness = nothing,
            local_unit_witnesses = (;),
            split_witness = nothing,
            bezout_witness = nothing,
        )

    recomputed = selection_ok ?
        _sl3_murthy_quillen_local_evidence_provider_fields(
            provider.witness_selection;
            witness = provider_inputs.witness,
            local_unit_witnesses = provider_inputs.local_unit_witnesses,
            split_witness = provider_inputs.split_witness,
            bezout_witness = provider_inputs.bezout_witness,
            metadata = provider.metadata,
        ) :
        nothing

    murthy_context_ok =
        recomputed !== nothing &&
        verify_sl3_local_murthy_input_context(provider.murthy_context) &&
        _same_sl3_local_murthy_input_context_data(
            provider.murthy_context,
            recomputed.murthy_context,
        )
    certificate_ok =
        recomputed !== nothing &&
        verify_sl3_local_realization(provider.murthy_certificate) &&
        _same_sl3_local_realization_certificate_data(
            provider.murthy_certificate,
            recomputed.murthy_certificate,
        )
    adapter_ok =
        recomputed !== nothing &&
        _verify_murthy_quillen_local_adapter(provider.murthy_adapter) &&
        _same_murthy_quillen_local_adapter_data(
            provider.murthy_adapter,
            recomputed.murthy_adapter,
        )
    selected_variable_ok =
        recomputed !== nothing &&
        provider.selected_variable == recomputed.selected_variable
    selected_variable_index_ok =
        recomputed !== nothing &&
        provider.selected_variable_index == recomputed.selected_variable_index
    local_product_ok =
        recomputed !== nothing &&
        provider.local_product == recomputed.local_product
    denominator_metadata_ok =
        recomputed !== nothing &&
        provider.denominator_metadata == recomputed.denominator_metadata
    witness_metadata_ok =
        recomputed !== nothing &&
        provider.witness_metadata == recomputed.witness_metadata
    replay_metadata_ok =
        recomputed !== nothing &&
        provider.replay_metadata == recomputed.replay_metadata
    staged_diagnostic_ok =
        recomputed !== nothing &&
        provider.staged_diagnostic == recomputed.staged_diagnostic
    metadata_ok =
        recomputed !== nothing &&
        provider.metadata == recomputed.metadata
    expected_quillen_local_sequences =
        adapter_ok &&
        provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence ?
        quillen_local_sequences_from_murthy_adapters(
            provider.local_product,
            provider.selected_variable,
            [provider.murthy_adapter],
        ) :
        Any[]
    quillen_local_sequences_ok =
        recomputed !== nothing &&
        provider.quillen_local_sequences isa AbstractVector &&
        (
            provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence ?
            all(
                verify_quillen_local_factor_sequence_certificate,
                provider.quillen_local_sequences,
            ) &&
            _same_quillen_local_factor_sequence_certificates(
                provider.quillen_local_sequences,
                expected_quillen_local_sequences,
            ) :
            isempty(provider.quillen_local_sequences)
        )

    overall_core_ok =
        selection_ok &&
        context_ok &&
        murthy_context_ok &&
        certificate_ok &&
        adapter_ok &&
        selected_variable_ok &&
        selected_variable_index_ok &&
        local_product_ok &&
        denominator_metadata_ok &&
        witness_metadata_ok &&
        replay_metadata_ok &&
        staged_diagnostic_ok &&
        metadata_ok &&
        quillen_local_sequences_ok

    return (;
        selection_ok,
        context_ok,
        murthy_context_ok,
        certificate_ok,
        adapter_ok,
        selected_variable_ok,
        selected_variable_index_ok,
        local_product_ok,
        denominator_metadata_ok,
        witness_metadata_ok,
        replay_metadata_ok,
        staged_diagnostic_ok,
        metadata_ok,
        quillen_local_sequences_ok,
        overall_core_ok,
    )
end

function _sl3_murthy_quillen_local_evidence_provider_verification(provider)
    core = _sl3_murthy_quillen_local_evidence_provider_core_verification(provider)
    stored_verification_ok = provider.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_sl3_murthy_quillen_local_evidence_provider(provider)::Bool
    try
        return _sl3_murthy_quillen_local_evidence_provider_verification(provider).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sl3_murthy_quillen_local_evidence_provider(
    selection::SL3LocalFormWitnessSelection;
    witness = nothing,
    local_unit_witnesses = (;),
    split_witness = nothing,
    bezout_witness = nothing,
    metadata = (;),
)
    fields = _sl3_murthy_quillen_local_evidence_provider_fields(
        selection;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
        metadata,
    )
    unchecked = SL3MurthyQuillenLocalEvidenceProvider(
        values(merge(fields, (; verification = nothing,)))...,
    )
    verification = _sl3_murthy_quillen_local_evidence_provider_core_verification(unchecked)
    checked = SL3MurthyQuillenLocalEvidenceProvider(
        values(merge(fields, (; verification,)))...,
    )
    _verify_sl3_murthy_quillen_local_evidence_provider(checked) ||
        error("internal SL_3 Murthy/Quillen local evidence provider verification failed")
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

function _polynomial_sl3_quillen_murthy_staged_failure_message(sl3_error)
    detail = sl3_error === nothing ? "" : ": $(sprint(showerror, sl3_error))"
    return "determinant-one polynomial input is outside the implemented evidence-backed SL_3 polynomial route$(detail)"
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

    if length(collect(gens(base_ring(A)))) > 1
        sl3_error = n == 3 ? _polynomial_sl3_quillen_murthy_route_error(A) : nothing
        throw(ArgumentError(_polynomial_sl3_quillen_murthy_staged_failure_message(sl3_error)))
    end

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

        if n > 3 && allow_recursive_column_peel
            recursive_staged_evidence =
                _polynomial_recursive_column_peel_public_staged_failure_evidence(A)
            if recursive_staged_evidence === nothing
                return _polynomial_recursive_column_peel_route_certificate(
                    A;
                    route_tag = :polynomial_column_peel,
                )
            elseif recursive_staged_evidence !== _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE
                return _polynomial_staged_failure_route_certificate(
                    A;
                    allow_recursive_column_peel = allow_recursive_column_peel,
                )
            end
        end

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

        if n == 3
            try
                return _polynomial_sl3_quillen_murthy_route_certificate(A)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
            end
        end

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
        return _polynomial_staged_failure_route_certificate(
            A;
            allow_recursive_column_peel = allow_recursive_column_peel,
        )
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
    evidence = merge(evidence, (; allow_recursive_column_peel))
    return _polynomial_route_certificate(A, :staged_failure, factors, product, evidence, :staged)
end

function _polynomial_recursive_column_peel_route_certificate(
    A;
    route_tag::Symbol = :polynomial_column_peel,
)
    _is_polynomial_column_peel_route(route_tag) ||
        throw(ArgumentError("unsupported polynomial column-peel route tag $(route_tag)"))
    evidence = _polynomial_column_peel_certificate(A)
    if evidence.final_certificate.evidence isa PolynomialSL3IdentityQuillenRouteEvidence
        throw(ArgumentError("polynomial column-peel route with a trivial identity SL_3 final block is certificate-only and not part of public elementary factorization dispatch"))
    end
    if !_polynomial_column_peel_public_mainline_supported(evidence)
        reason_code = _polynomial_column_peel_public_reason_code(evidence)
        throw(ArgumentError(_polynomial_column_peel_public_staged_message(reason_code)))
    end
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

function _polynomial_sl3_supplied_quillen_route_fields(A; metadata = (;))
    nrows(A) == 3 && ncols(A) == 3 ||
        throw(ArgumentError("SL_3 supplied-Quillen route requires a 3 x 3 input"))
    data = _polynomial_quillen_supplied_evidence_data(A)
    data === nothing &&
        throw(ArgumentError("SL_3 supplied-Quillen route requires supplied Quillen evidence"))
    patch = _polynomial_quillen_supplied_evidence_patch(A)
    adapter = _polynomial_quillen_patch_route_adapter(A, patch)
    replay_metadata = (;
        source = :sl3_supplied_quillen_polynomial_route,
        route_issue_id = "#263",
        driver_issue_id = "#184",
        patch_replay_metadata = adapter.quillen_patch.replay_metadata,
        adapter_replay_metadata = adapter.replay_metadata,
        metadata,
    )
    return (;
        target = adapter.target_matrix,
        route = :quillen_patch,
        supplied_evidence = data,
        quillen_route_adapter = adapter,
        replay_metadata,
    )
end

function _same_polynomial_sl3_supplied_quillen_data(left, right)::Bool
    local_certificates_ok =
        _same_quillen_local_factor_sequence_certificates(
            left.local_certificates,
            right.local_certificates,
        ) ||
        _same_polynomial_sl3_supplied_quillen_local_certificates(
            left.local_certificates,
            right.local_certificates,
        )
    return left.selected_variable == right.selected_variable &&
        left.cover_generator == right.cover_generator &&
        left.row == right.row &&
        left.col == right.col &&
        left.entry == right.entry &&
        left.base_entry == right.base_entry &&
        left.delta_entry == right.delta_entry &&
        left.base_term_policy == right.base_term_policy &&
        _polynomial_route_factor_sequences_equal(left.base_term_factors, right.base_term_factors) &&
        local_certificates_ok &&
        all(verify_quillen_local_factor_sequence_certificate, left.local_certificates) &&
        all(verify_quillen_local_factor_sequence_certificate, right.local_certificates)
end

function _same_polynomial_sl3_supplied_quillen_elementary_correction(left, right)::Bool
    return left.row == right.row &&
        left.col == right.col &&
        left.entry == right.entry
end

function _same_polynomial_sl3_supplied_quillen_local_contribution(left, right)::Bool
    return left.certificate.indices == right.certificate.indices &&
        left.certificate.denominators == right.certificate.denominators &&
        left.denominator == right.denominator &&
        left.coverage_multiplier == right.coverage_multiplier &&
        _same_polynomial_sl3_supplied_quillen_elementary_correction(
            left.correction,
            right.correction,
        )
end

function _same_polynomial_sl3_supplied_quillen_local_contributions(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _same_polynomial_sl3_supplied_quillen_local_contribution(
            left[idx],
            right[idx],
        ) || return false
    end
    return true
end

function _same_polynomial_sl3_supplied_quillen_local_certificate(left, right)::Bool
    return left.original_input == right.original_input &&
        left.ring == right.ring &&
        left.size == right.size &&
        left.selected_variable == right.selected_variable &&
        _same_quillen_local_elementary_factors(left.factors, right.factors) &&
        left.raw_denominators == right.raw_denominators &&
        left.product_denominator == right.product_denominator &&
        left.local_product == right.local_product &&
        left.local_correction == right.local_correction &&
        _same_polynomial_sl3_supplied_quillen_local_contributions(
            left.normalized_local_contributions,
            right.normalized_local_contributions,
        ) &&
        _same_quillen_factors(
            left.normalized_global_elementary_factors,
            right.normalized_global_elementary_factors,
        ) &&
        left.patched_substitution_witness == right.patched_substitution_witness &&
        left.chain_witness == right.chain_witness &&
        left.witness_metadata == right.witness_metadata
end

function _same_polynomial_sl3_supplied_quillen_local_certificates(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _same_polynomial_sl3_supplied_quillen_local_certificate(
            left[idx],
            right[idx],
        ) || return false
    end
    return true
end

function _same_polynomial_sl3_supplied_quillen_cover_candidate_verification(left, right)::Bool
    return left.local_count == right.local_count &&
        left.raw_denominators == right.raw_denominators &&
        left.local_certificates_ok == right.local_certificates_ok &&
        left.same_original_input_ok == right.same_original_input_ok &&
        left.same_ring_ok == right.same_ring_ok &&
        left.same_size_ok == right.same_size_ok &&
        left.same_selected_variable_ok == right.same_selected_variable_ok &&
        _same_quillen_local_denominator_supports(left.local_supports, right.local_supports) &&
        left.local_supports_ok == right.local_supports_ok &&
        left.raw_denominators_ok == right.raw_denominators_ok &&
        left.replay_metadata_ok == right.replay_metadata_ok &&
        left.overall_ok == right.overall_ok
end

function _same_polynomial_sl3_supplied_quillen_cover_candidate(left, right)::Bool
    return _same_quillen_denominator_cover_candidate_data(left, right) ||
        (
            left.original_input == right.original_input &&
            left.ring == right.ring &&
            left.size == right.size &&
            left.selected_variable == right.selected_variable &&
            _same_polynomial_sl3_supplied_quillen_local_certificates(
                left.local_certificates,
                right.local_certificates,
            ) &&
            left.raw_denominators == right.raw_denominators &&
            _same_quillen_local_denominator_supports(
                left.local_supports,
                right.local_supports,
            ) &&
            _same_polynomial_sl3_supplied_quillen_cover_candidate_verification(
                left.verification,
                right.verification,
            )
        )
end

function _same_polynomial_sl3_supplied_quillen_solver_result(left, right)::Bool
    return _same_quillen_denominator_cover_solver_result_data(left, right) ||
        (
            _same_polynomial_sl3_supplied_quillen_cover_candidate(
                left.source_candidate,
                right.source_candidate,
            ) &&
            left.ring == right.ring &&
            left.raw_denominators == right.raw_denominators &&
            left.exponent == right.exponent &&
            left.powered_denominators == right.powered_denominators &&
            left.coverage_multipliers == right.coverage_multipliers &&
            left.coverage_terms == right.coverage_terms &&
            left.coverage_sum == right.coverage_sum &&
            _same_quillen_cover_certificate_data(
                left.cover_certificate,
                right.cover_certificate,
            ) &&
            _same_quillen_denominator_cover_solver_verification(
                left.verification,
                right.verification,
            )
        )
end

function _same_polynomial_sl3_supplied_quillen_sequence_expansion_verification(left, right)::Bool
    return left.local_certificate_ok == right.local_certificate_ok &&
        left.solver_result_ok == right.solver_result_ok &&
        left.local_index_ok == right.local_index_ok &&
        left.solver_context_ok == right.solver_context_ok &&
        left.powered_denominator == right.powered_denominator &&
        left.coverage_multiplier == right.coverage_multiplier &&
        left.cover_term == right.cover_term &&
        left.cover_term_ok == right.cover_term_ok &&
        left.factor_provenance == right.factor_provenance &&
        left.factor_provenance_ok == right.factor_provenance_ok &&
        _same_quillen_factors(left.global_elementary_factors, right.global_elementary_factors) &&
        left.global_elementary_factors_ok == right.global_elementary_factors_ok &&
        left.replay_metadata_ok == right.replay_metadata_ok &&
        left.overall_ok == right.overall_ok
end

function _same_polynomial_sl3_supplied_quillen_sequence_expansions(left, right)::Bool
    return _same_quillen_sequence_expansions(left, right) ||
        (
            length(left) == length(right) &&
            all(eachindex(left, right)) do idx
                left[idx].local_index == right[idx].local_index &&
                _same_polynomial_sl3_supplied_quillen_solver_result(
                    left[idx].solver_result,
                    right[idx].solver_result,
                ) &&
                left[idx].cover_term == right[idx].cover_term &&
                _same_quillen_factors(
                    left[idx].global_elementary_factors,
                    right[idx].global_elementary_factors,
                ) &&
                _same_polynomial_sl3_supplied_quillen_sequence_expansion_verification(
                    left[idx].verification,
                    right[idx].verification,
                )
            end
        )
end

function _same_polynomial_sl3_supplied_quillen_substitution_chain_verification(left, right)::Bool
    return left.solver_result_ok == right.solver_result_ok &&
        left.ring_ok == right.ring_ok &&
        left.matrix_ok == right.matrix_ok &&
        left.selected_variable_ok == right.selected_variable_ok &&
        left.sign_convention_ok == right.sign_convention_ok &&
        left.coefficient_count == right.coefficient_count &&
        left.step_count == right.step_count &&
        left.bracket_count == right.bracket_count &&
        left.cumulative_coefficients == right.cumulative_coefficients &&
        left.cumulative_coefficients_ok == right.cumulative_coefficients_ok &&
        left.intermediate_matrices == right.intermediate_matrices &&
        left.intermediate_matrices_ok == right.intermediate_matrices_ok &&
        _same_quillen_patch_substitution_steps(left.expected_steps, right.expected_steps) &&
        left.steps_ok == right.steps_ok &&
        left.bracket_matrices == right.bracket_matrices &&
        left.bracket_matrices_ok == right.bracket_matrices_ok &&
        left.final_coefficient == right.final_coefficient &&
        left.final_coefficient_ok == right.final_coefficient_ok &&
        left.base_term == right.base_term &&
        left.base_term_ok == right.base_term_ok &&
        left.telescope_product == right.telescope_product &&
        left.telescope_ok == right.telescope_ok &&
        left.replay_metadata_ok == right.replay_metadata_ok &&
        left.overall_ok == right.overall_ok
end

function _same_polynomial_sl3_supplied_quillen_substitution_chain(left, right)::Bool
    return _same_polynomial_sl3_quillen_murthy_substitution_chain_data(left, right) ||
        (
            left.original_matrix == right.original_matrix &&
            left.ring == right.ring &&
            left.size == right.size &&
            left.selected_variable == right.selected_variable &&
            left.sign_convention == right.sign_convention &&
            _same_polynomial_sl3_supplied_quillen_solver_result(
                left.solver_result,
                right.solver_result,
            ) &&
            left.cumulative_coefficients == right.cumulative_coefficients &&
            left.intermediate_matrices == right.intermediate_matrices &&
            _same_quillen_patch_substitution_steps(left.steps, right.steps) &&
            left.bracket_matrices == right.bracket_matrices &&
            left.base_term == right.base_term &&
            left.metadata == right.metadata &&
            _same_polynomial_sl3_supplied_quillen_substitution_chain_verification(
                left.verification,
                right.verification,
            )
        )
end

function _same_polynomial_sl3_supplied_quillen_patch_verification(left, right)::Bool
    return _same_quillen_supplied_patch_verification(left, right) ||
        (
            left.local_certificates_ok == right.local_certificates_ok &&
            left.denominator_candidate_ok == right.denominator_candidate_ok &&
            left.denominator_candidate_matches == right.denominator_candidate_matches &&
            left.solver_result_ok == right.solver_result_ok &&
            left.solver_source_candidate_ok == right.solver_source_candidate_ok &&
            left.cover_certificate_ok == right.cover_certificate_ok &&
            left.substitution_chain_ok == right.substitution_chain_ok &&
            left.substitution_chain_matches == right.substitution_chain_matches &&
            left.base_term_ok == right.base_term_ok &&
            left.sequence_expansions_ok == right.sequence_expansions_ok &&
            _same_quillen_factors(
                left.global_elementary_factors,
                right.global_elementary_factors,
            ) &&
            left.global_elementary_factors_ok == right.global_elementary_factors_ok &&
            left.product == right.product &&
            left.product_ok == right.product_ok &&
            left.target == right.target &&
            left.target_ok == right.target_ok &&
            left.replay_metadata_ok == right.replay_metadata_ok &&
            left.overall_ok == right.overall_ok
        )
end

function _same_polynomial_sl3_supplied_quillen_patch_data(adapted, raw)::Bool
    return verify_quillen_patch(adapted) &&
        verify_quillen_patch(raw) &&
        adapted.ring == raw.ring &&
        adapted.size == raw.size &&
        adapted.substitution_variable == raw.substitution_variable &&
        adapted.original_input == raw.original_input &&
        (
            _same_quillen_local_factor_sequence_certificates(
                adapted.local_certificates,
                raw.local_certificates,
            ) ||
            _same_polynomial_sl3_supplied_quillen_local_certificates(
                adapted.local_certificates,
                raw.local_certificates,
            )
        ) &&
        _same_polynomial_sl3_supplied_quillen_cover_candidate(
            adapted.denominator_candidate,
            raw.denominator_candidate,
        ) &&
        _same_polynomial_sl3_supplied_quillen_solver_result(
            adapted.solver_result,
            raw.solver_result,
        ) &&
        _same_quillen_cover_certificate_data(
            adapted.cover_certificate,
            raw.cover_certificate,
        ) &&
        _same_polynomial_sl3_supplied_quillen_substitution_chain(
            adapted.substitution_chain,
            raw.substitution_chain,
        ) &&
        adapted.base_term_policy == raw.base_term_policy &&
        adapted.base_term == raw.base_term &&
        _polynomial_route_factor_sequences_equal(
            adapted.base_term_factors,
            raw.base_term_factors,
        ) &&
        adapted.base_term_product == raw.base_term_product &&
        _same_polynomial_sl3_supplied_quillen_sequence_expansions(
            adapted.sequence_expansions,
            raw.sequence_expansions,
        ) &&
        _same_quillen_factors(
            adapted.sequence_elementary_factors,
            raw.sequence_elementary_factors,
        ) &&
        _same_quillen_factors(
            adapted.global_elementary_factors,
            raw.global_elementary_factors,
        ) &&
        adapted.product == raw.product &&
        adapted.target == raw.target &&
        _same_polynomial_sl3_supplied_quillen_patch_verification(
            adapted.verification,
            raw.verification,
        )
end

function _same_polynomial_sl3_supplied_quillen_route_adapter(left, right)::Bool
    return left.target == right.target &&
        left.route == right.route &&
        _same_polynomial_sl3_supplied_quillen_patch_data(
            left.quillen_patch,
            right.quillen_patch,
        ) &&
        left.target_matrix == right.target_matrix &&
        _polynomial_route_factor_sequences_equal(
            left.global_elementary_factors,
            right.global_elementary_factors,
        ) &&
        left.product == right.product &&
        left.replay_metadata.patch_size == right.replay_metadata.patch_size &&
        left.replay_metadata.substitution_variable == right.replay_metadata.substitution_variable &&
        left.replay_metadata.denominator_data == right.replay_metadata.denominator_data &&
        left.replay_metadata.local_certificate_count == right.replay_metadata.local_certificate_count &&
        left.replay_metadata.normalized_contribution_count ==
            right.replay_metadata.normalized_contribution_count &&
        _verify_polynomial_quillen_patch_route_adapter(left) &&
        _verify_polynomial_quillen_patch_route_adapter(right)
end

function _polynomial_sl3_supplied_quillen_route_core_verification(evidence)
    try
        route_ok = evidence.route == :quillen_patch
        fields = route_ok ? _polynomial_sl3_supplied_quillen_route_fields(evidence.target) : nothing
        supplied_evidence_ok =
            fields !== nothing &&
            _same_polynomial_sl3_supplied_quillen_data(
                evidence.supplied_evidence,
                fields.supplied_evidence,
            )
        adapter_ok =
            fields !== nothing &&
            _same_polynomial_sl3_supplied_quillen_route_adapter(
                evidence.quillen_route_adapter,
                fields.quillen_route_adapter,
            )
        replay_metadata_ok =
            evidence.replay_metadata.source == :sl3_supplied_quillen_polynomial_route &&
            evidence.replay_metadata.route_issue_id == "#263" &&
            evidence.replay_metadata.driver_issue_id == "#184" &&
            evidence.replay_metadata.patch_replay_metadata ==
                evidence.quillen_route_adapter.quillen_patch.replay_metadata &&
            evidence.replay_metadata.adapter_replay_metadata ==
                evidence.quillen_route_adapter.replay_metadata
        overall_core_ok = route_ok && supplied_evidence_ok && adapter_ok && replay_metadata_ok
        return (;
            route_ok,
            supplied_evidence_ok,
            adapter_ok,
            replay_metadata_ok,
            overall_core_ok,
        )
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return (; # COV_EXCL_LINE
            route_ok = false,
            supplied_evidence_ok = false,
            adapter_ok = false,
            replay_metadata_ok = false,
            overall_core_ok = false,
        )
    end
end

function _polynomial_sl3_supplied_quillen_route_verification(evidence)
    core = _polynomial_sl3_supplied_quillen_route_core_verification(evidence)
    stored_verification_ok = evidence.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_polynomial_sl3_supplied_quillen_route_evidence(evidence)::Bool
    try
        return _polynomial_sl3_supplied_quillen_route_verification(evidence).overall_ok
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return false # COV_EXCL_LINE
    end
end

function _polynomial_sl3_supplied_quillen_route_evidence(A; metadata = (;))
    fields = _polynomial_sl3_supplied_quillen_route_fields(A; metadata)
    raw = PolynomialSL3SuppliedQuillenRouteEvidence(
        values(merge(fields, (; verification = nothing,)))...,
    )
    verification = _polynomial_sl3_supplied_quillen_route_core_verification(raw)
    evidence = PolynomialSL3SuppliedQuillenRouteEvidence(
        values(merge(fields, (; verification,)))...,
    )
    _verify_polynomial_sl3_supplied_quillen_route_evidence(evidence) ||
        error("internal SL_3 supplied-Quillen route evidence verification failed")
    return evidence
end

function _polynomial_sl3_supplied_quillen_route_certificate(A)
    evidence = _polynomial_sl3_supplied_quillen_route_evidence(A)
    adapter = evidence.quillen_route_adapter
    factors = copy(adapter.global_elementary_factors)
    return _polynomial_route_certificate(
        adapter.target_matrix,
        :quillen_patch,
        factors,
        adapter.product,
        evidence,
        :supported,
    )
end

function _polynomial_sl3_identity_quillen_route_metadata(A)
    return (;
        source = :identity_sl3_quillen_route,
        route = :quillen_patch,
        issue_refs = ("#184", "#263"),
        size = nrows(A),
        factor_count = 0,
        proof = :empty_elementary_factorization,
    )
end

function _polynomial_sl3_identity_quillen_route_core_verification(evidence)
    try
        target = evidence.target
        route_ok = evidence.route == :quillen_patch
        square_ok = nrows(target) == 3 && ncols(target) == 3
        R = base_ring(target)
        ring_profile_ok = _factorization_ring_profile(R) == :polynomial
        coefficient_ring_ok = _polynomial_exact_field_backed_ring(R)
        identity_target_ok = target == identity_matrix(R, 3)
        factors_vector_ok = evidence.factors isa AbstractVector
        factors_empty_ok = factors_vector_ok && isempty(evidence.factors)
        replayed_product =
            factors_vector_ok ? _polynomial_route_factor_product(evidence.factors, R, 3) : nothing
        product_replay_ok = replayed_product !== nothing && replayed_product == identity_matrix(R, 3)
        product_ok = evidence.product == identity_matrix(R, 3)
        replay_metadata_ok =
            evidence.replay_metadata == _polynomial_sl3_identity_quillen_route_metadata(target)
        factorization_ok = factors_vector_ok && verify_factorization(target, evidence.factors)
        overall_core_ok =
            route_ok &&
            square_ok &&
            ring_profile_ok &&
            coefficient_ring_ok &&
            identity_target_ok &&
            factors_empty_ok &&
            product_replay_ok &&
            product_ok &&
            replay_metadata_ok &&
            factorization_ok
        return (;
            route_ok,
            square_ok,
            ring_profile_ok,
            coefficient_ring_ok,
            identity_target_ok,
            factors_vector_ok,
            factors_empty_ok,
            product_replay_ok,
            product_ok,
            replay_metadata_ok,
            factorization_ok,
            overall_core_ok,
        )
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return (; # COV_EXCL_LINE
            route_ok = false,
            square_ok = false,
            ring_profile_ok = false,
            coefficient_ring_ok = false,
            identity_target_ok = false,
            factors_vector_ok = false,
            factors_empty_ok = false,
            product_replay_ok = false,
            product_ok = false,
            replay_metadata_ok = false,
            factorization_ok = false,
            overall_core_ok = false,
        )
    end
end

function _polynomial_sl3_identity_quillen_route_verification(evidence)
    core = _polynomial_sl3_identity_quillen_route_core_verification(evidence)
    stored_verification_ok = evidence.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_polynomial_sl3_identity_quillen_route_evidence(evidence)::Bool
    try
        return _polynomial_sl3_identity_quillen_route_verification(evidence).overall_ok
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return false # COV_EXCL_LINE
    end
end

function _polynomial_sl3_identity_quillen_route_evidence(A)
    nrows(A) == 3 && ncols(A) == 3 ||
        throw(ArgumentError("SL_3 identity Quillen route requires a 3 x 3 input"))
    R = base_ring(A)
    _factorization_ring_profile(R) == :polynomial ||
        throw(ArgumentError("SL_3 identity Quillen route requires an ordinary polynomial ring"))
    _polynomial_exact_field_backed_ring(R) ||
        throw(ArgumentError("SL_3 identity Quillen route requires an exact field-backed coefficient ring"))
    A == identity_matrix(R, 3) ||
        throw(ArgumentError("SL_3 identity Quillen route requires the identity matrix"))

    product = identity_matrix(R, 3)
    factors = typeof(product)[]
    replay_metadata = _polynomial_sl3_identity_quillen_route_metadata(A)
    raw = PolynomialSL3IdentityQuillenRouteEvidence(
        A,
        :quillen_patch,
        factors,
        product,
        replay_metadata,
        nothing,
    )
    verification = _polynomial_sl3_identity_quillen_route_core_verification(raw)
    evidence = PolynomialSL3IdentityQuillenRouteEvidence(
        A,
        :quillen_patch,
        factors,
        product,
        replay_metadata,
        verification,
    )
    _verify_polynomial_sl3_identity_quillen_route_evidence(evidence) ||
        error("internal SL_3 identity Quillen route evidence verification failed")
    return evidence
end

function _polynomial_sl3_identity_quillen_route_certificate(A)
    evidence = _polynomial_sl3_identity_quillen_route_evidence(A)
    factors = copy(evidence.factors)
    return _polynomial_route_certificate(
        evidence.target,
        :quillen_patch,
        factors,
        evidence.product,
        evidence,
        :supported,
    )
end

function _polynomial_sl3_quillen_murthy_selected_variable(A)
    R = base_ring(A)
    entries = _sl3_local_target_entries(A)
    entries === nothing &&
        throw(ArgumentError("SL_3 Quillen/Murthy route is missing a #236 local-form witness: input is not an already-special-form SL_3 context"))
    for (index, generator) in enumerate(gens(R))
        monicity_witness = _sl3_local_monicity_witness(entries.p, index, R)
        if monicity_witness.is_monic && monicity_witness.degree > 0
            return generator, index, entries
        end
    end
    throw(ArgumentError("SL_3 Quillen/Murthy route is missing a #236 local-form witness: special-form p is not monic in any selected generator"))
end

function _polynomial_sl3_quillen_murthy_local_form_witness(A, selected_variable, entries)
    return (;
        entries,
        p = entries.p,
        monic_entry_position = (1, 1),
        source_matrix = A,
        selected_variable,
        replay_steps = ((; kind = :issue238_already_special_form_replay),),
        witness_issue_id = "#236",
        route_issue_id = "#238",
    )
end

function _polynomial_sl3_quillen_murthy_route_metadata(
    context,
    selection,
    provider,
    base_term_policy::Symbol,
    base_term_factors,
    metadata,
)
    return (;
        source = :sl3_quillen_murthy_polynomial_route,
        route_issue_id = "#238",
        context_issue_id = "#235",
        witness_issue_id = "#236",
        provider_issue_id = "#237",
        adapter_issue_id = "#219",
        patch_issue_id = "#220",
        context_metadata = context.catalog_metadata,
        witness_metadata = selection.local_form_witness,
        provider_replay_metadata = provider.replay_metadata,
        base_term_policy,
        base_term_factor_count = length(base_term_factors),
        metadata,
    )
end

function _polynomial_sl3_quillen_murthy_route_patch(
    patch,
    route_metadata,
)
    replay_metadata = _quillen_supplied_patch_metadata(
        patch.denominator_candidate,
        patch.solver_result,
        patch.substitution_chain,
        patch.base_term_policy,
        patch.sequence_expansions,
        route_metadata,
    )
    provisional = QuillenSuppliedEvidencePatchAssembly(
        patch.ring,
        patch.size,
        patch.substitution_variable,
        patch.original_input,
        patch.local_certificates,
        patch.denominator_candidate,
        patch.solver_result,
        patch.cover_certificate,
        patch.substitution_chain,
        patch.base_term_policy,
        patch.base_term,
        patch.base_term_factors,
        patch.base_term_product,
        patch.sequence_expansions,
        patch.sequence_elementary_factors,
        patch.global_elementary_factors,
        patch.product,
        patch.target,
        replay_metadata,
        patch.verification,
    )
    verification = replay_quillen_supplied_evidence_patch(provisional)
    checked = QuillenSuppliedEvidencePatchAssembly(
        patch.ring,
        patch.size,
        patch.substitution_variable,
        patch.original_input,
        patch.local_certificates,
        patch.denominator_candidate,
        patch.solver_result,
        patch.cover_certificate,
        patch.substitution_chain,
        patch.base_term_policy,
        patch.base_term,
        patch.base_term_factors,
        patch.base_term_product,
        patch.sequence_expansions,
        patch.sequence_elementary_factors,
        patch.global_elementary_factors,
        patch.product,
        patch.target,
        replay_metadata,
        verification,
    )
    verify_quillen_patch(checked) ||
        error("internal SL_3 Quillen/Murthy route patch adaptation failed")
    return checked
end

function _polynomial_sl3_quillen_murthy_route_consumption_metadata(
    raw_consumption,
    route_metadata,
    patch,
)
    return (;
        source = :sl3_quillen_murthy_polynomial_route,
        route_issue_id = "#238",
        route_metadata,
        raw_consumption_replay_metadata = raw_consumption.replay_metadata,
        patch_replay_metadata = patch.replay_metadata,
    )
end

function _polynomial_sl3_quillen_murthy_raw_route_metadata(raw_consumption)
    hasproperty(raw_consumption.replay_metadata, :metadata) || return nothing
    return raw_consumption.replay_metadata.metadata
end

function _polynomial_sl3_quillen_murthy_route_metadata_payload(route_metadata)
    route_metadata === nothing && return (;)
    hasproperty(route_metadata, :metadata) || return (;)
    return route_metadata.metadata
end

function _same_polynomial_sl3_quillen_murthy_substitution_chain_data(left, right)::Bool
    return left.original_matrix == right.original_matrix &&
           left.ring == right.ring &&
           left.size == right.size &&
           left.selected_variable == right.selected_variable &&
           left.sign_convention == right.sign_convention &&
           _same_quillen_denominator_cover_solver_result_data(
               left.solver_result,
               right.solver_result,
           ) &&
           left.cumulative_coefficients == right.cumulative_coefficients &&
           left.intermediate_matrices == right.intermediate_matrices &&
           _same_quillen_patch_substitution_steps(left.steps, right.steps) &&
           left.bracket_matrices == right.bracket_matrices &&
           left.base_term == right.base_term &&
           left.metadata == right.metadata &&
           left.replay_metadata == right.replay_metadata &&
           _same_quillen_patch_substitution_chain_verification(
               left.verification,
               right.verification,
           )
end

function _same_polynomial_sl3_quillen_murthy_route_patch_verification(left, right)::Bool
    return left.local_certificates_ok == right.local_certificates_ok &&
           left.denominator_candidate_ok == right.denominator_candidate_ok &&
           left.denominator_candidate_matches == right.denominator_candidate_matches &&
           left.solver_result_ok == right.solver_result_ok &&
           left.solver_source_candidate_ok == right.solver_source_candidate_ok &&
           left.cover_certificate_ok == right.cover_certificate_ok &&
           left.substitution_chain_ok == right.substitution_chain_ok &&
           left.substitution_chain_matches == right.substitution_chain_matches &&
           left.base_term_ok == right.base_term_ok &&
           left.sequence_expansions_ok == right.sequence_expansions_ok &&
           _same_quillen_factors(
               left.global_elementary_factors,
               right.global_elementary_factors,
           ) &&
           left.global_elementary_factors_ok == right.global_elementary_factors_ok &&
           left.product == right.product &&
           left.product_ok == right.product_ok &&
           left.target == right.target &&
           left.target_ok == right.target_ok &&
           left.replay_metadata_ok &&
           right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _same_polynomial_sl3_quillen_murthy_rewritten_patch_data(adapted, raw)::Bool
    return verify_quillen_patch(adapted) &&
           verify_quillen_patch(raw) &&
           adapted.ring == raw.ring &&
           adapted.size == raw.size &&
           adapted.substitution_variable == raw.substitution_variable &&
           adapted.original_input == raw.original_input &&
           adapted.local_certificates == raw.local_certificates &&
           _same_quillen_local_factor_sequence_certificates(
               adapted.local_certificates,
               raw.local_certificates,
           ) &&
           _same_quillen_denominator_cover_candidate_data(
               adapted.denominator_candidate,
               raw.denominator_candidate,
           ) &&
           _same_quillen_denominator_cover_solver_result_data(
               adapted.solver_result,
               raw.solver_result,
           ) &&
           _same_quillen_cover_certificate_data(
               adapted.cover_certificate,
               raw.cover_certificate,
           ) &&
           _same_polynomial_sl3_quillen_murthy_substitution_chain_data(
               adapted.substitution_chain,
               raw.substitution_chain,
           ) &&
           adapted.base_term_policy == raw.base_term_policy &&
           adapted.base_term == raw.base_term &&
           _polynomial_route_factor_sequences_equal(
               adapted.base_term_factors,
               raw.base_term_factors,
           ) &&
           adapted.base_term_product == raw.base_term_product &&
           _same_quillen_sequence_expansions(
               adapted.sequence_expansions,
               raw.sequence_expansions,
           ) &&
           _same_quillen_factors(
               adapted.sequence_elementary_factors,
               raw.sequence_elementary_factors,
           ) &&
           _same_quillen_factors(
               adapted.global_elementary_factors,
               raw.global_elementary_factors,
           ) &&
           adapted.product == raw.product &&
           adapted.target == raw.target &&
           _same_polynomial_sl3_quillen_murthy_route_patch_verification(
               adapted.verification,
               raw.verification,
           )
end

function _polynomial_sl3_quillen_murthy_route_consumption_core_verification(consumption)
    raw_consumption_ok = verify_quillen_murthy_adapter_consumption(consumption.raw_consumption)
    adapter_context_ok =
        raw_consumption_ok &&
        consumption.original_input == consumption.raw_consumption.original_input &&
        consumption.ring == consumption.raw_consumption.ring &&
        consumption.size == consumption.raw_consumption.size &&
        consumption.selected_variable == consumption.raw_consumption.selected_variable
    adapters_ok =
        adapter_context_ok &&
        length(consumption.murthy_adapters) == length(consumption.raw_consumption.murthy_adapters) &&
        all(
            idx -> _same_murthy_quillen_local_adapter_data(
                consumption.murthy_adapters[idx],
                consumption.raw_consumption.murthy_adapters[idx],
            ),
            eachindex(consumption.murthy_adapters, consumption.raw_consumption.murthy_adapters),
        )
    local_sequences_ok =
        adapters_ok &&
        _same_quillen_local_factor_sequence_certificates(
            consumption.local_sequence_certificates,
            consumption.raw_consumption.local_sequence_certificates,
        ) &&
        all(
            verify_quillen_local_factor_sequence_certificate,
            consumption.local_sequence_certificates,
        )
    raw_route_metadata =
        raw_consumption_ok ?
        _polynomial_sl3_quillen_murthy_raw_route_metadata(consumption.raw_consumption) :
        nothing
    patch_rewrite_ok =
        local_sequences_ok &&
        _same_polynomial_sl3_quillen_murthy_rewritten_patch_data(
            consumption.patch,
            consumption.raw_consumption.patch,
        )
    patch_route_metadata_ok =
        patch_rewrite_ok &&
        raw_route_metadata !== nothing &&
        hasproperty(consumption.patch.replay_metadata, :metadata) &&
        consumption.patch.replay_metadata.metadata == raw_route_metadata
    patch_ok =
        patch_route_metadata_ok &&
        consumption.patch.original_input == consumption.original_input &&
        consumption.patch.ring == consumption.ring &&
        consumption.patch.size == consumption.size &&
        consumption.patch.substitution_variable == consumption.selected_variable &&
        consumption.patch.base_term_policy == :already_handled &&
        _same_quillen_local_factor_sequence_certificates(
            consumption.patch.local_certificates,
            consumption.local_sequence_certificates,
        ) &&
        _polynomial_route_factor_sequences_equal(
            consumption.patch.global_elementary_factors,
            consumption.raw_consumption.patch.global_elementary_factors,
        )
    expected_replay_metadata =
        patch_ok ?
        _polynomial_sl3_quillen_murthy_route_consumption_metadata(
            consumption.raw_consumption,
            raw_route_metadata,
            consumption.patch,
        ) :
        nothing
    replay_metadata_ok = consumption.replay_metadata == expected_replay_metadata
    overall_core_ok =
        raw_consumption_ok &&
        adapter_context_ok &&
        adapters_ok &&
        local_sequences_ok &&
        patch_ok &&
        replay_metadata_ok
    return (;
        raw_consumption_ok,
        adapter_context_ok,
        adapters_ok,
        local_sequences_ok,
        patch_rewrite_ok,
        patch_route_metadata_ok,
        patch_ok,
        replay_metadata_ok,
        overall_core_ok,
    )
end

function _polynomial_sl3_quillen_murthy_route_consumption_verification(consumption)
    core = _polynomial_sl3_quillen_murthy_route_consumption_core_verification(consumption)
    stored_verification_ok = consumption.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function verify_quillen_murthy_adapter_consumption(
    consumption::PolynomialSL3QuillenMurthyRouteConsumption,
)::Bool
    try
        return _polynomial_sl3_quillen_murthy_route_consumption_verification(consumption).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_sl3_quillen_murthy_route_consumption(
    raw_consumption,
    route_metadata,
)
    patch = _polynomial_sl3_quillen_murthy_route_patch(raw_consumption.patch, route_metadata)
    replay_metadata = _polynomial_sl3_quillen_murthy_route_consumption_metadata(
        raw_consumption,
        route_metadata,
        patch,
    )
    fields = (;
        original_input = raw_consumption.original_input,
        ring = raw_consumption.ring,
        size = raw_consumption.size,
        selected_variable = raw_consumption.selected_variable,
        murthy_adapters = raw_consumption.murthy_adapters,
        local_sequence_certificates = raw_consumption.local_sequence_certificates,
        patch,
        raw_consumption,
        replay_metadata,
    )
    unchecked = PolynomialSL3QuillenMurthyRouteConsumption(
        values(merge(fields, (; verification = nothing,)))...,
    )
    verification = _polynomial_sl3_quillen_murthy_route_consumption_core_verification(unchecked)
    checked = PolynomialSL3QuillenMurthyRouteConsumption(
        values(merge(fields, (; verification,)))...,
    )
    verify_quillen_murthy_adapter_consumption(checked) ||
        error("internal SL_3 Quillen/Murthy route consumption verification failed")
    return checked
end

function _polynomial_sl3_quillen_murthy_route_fields(A; metadata = (;))
    nrows(A) == 3 && ncols(A) == 3 ||
        throw(ArgumentError("SL_3 Quillen/Murthy route requires a 3 x 3 input"))
    R = base_ring(A)
    _factorization_ring_profile(R) == :polynomial ||
        throw(ArgumentError("SL_3 Quillen/Murthy route requires an ordinary polynomial ring"))
    _polynomial_exact_field_backed_ring(R) ||
        throw(ArgumentError(_polynomial_unsupported_coefficient_ring_message()))
    _require_polynomial_sl_determinant(A)
    selected_variable, selected_variable_index, entries =
        _polynomial_sl3_quillen_murthy_selected_variable(A)
    local_form_witness =
        _polynomial_sl3_quillen_murthy_local_form_witness(A, selected_variable, entries)
    context = _sl3_realization_input_context(
        A;
        selected_variable = (;
            name = string(selected_variable),
            generator = selected_variable,
            index = selected_variable_index,
            status = :passes,
        ),
        catalog_metadata = (;
            source = :sl3_quillen_murthy_polynomial_route,
            route_issue_id = "#238",
            context_issue_id = "#235",
            driver_issue_id = "#184",
        ),
        local_form_witness,
    )
    selection = _select_sl3_local_form_witness(context)
    provider = try
        _sl3_murthy_quillen_local_evidence_provider(
            selection;
            metadata = (;
                source = :sl3_quillen_murthy_polynomial_route,
                route_issue_id = "#238",
                provider_issue_id = "#237",
            ),
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        throw(ArgumentError("SL_3 Quillen/Murthy route is missing #237 ordinary Quillen local evidence: $(sprint(showerror, err))"))
    end
    provider.staged_diagnostic.status == :supported ||
        throw(ArgumentError("SL_3 Quillen/Murthy route is missing #237 ordinary Quillen local evidence: $(provider.staged_diagnostic.message)"))
    provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence ||
        throw(ArgumentError("SL_3 Quillen/Murthy route requires ordinary Quillen local sequence evidence before #220 patching"))
    base_term_factors = typeof(A)[]
    base_term_policy = :already_handled
    route_metadata = _polynomial_sl3_quillen_murthy_route_metadata(
        context,
        selection,
        provider,
        base_term_policy,
        base_term_factors,
        metadata,
    )
    raw_consumption = try
        consume_murthy_quillen_adapters_for_patch(
            A,
            selected_variable,
            [provider.murthy_adapter];
            base_term_policy,
            base_term_factors,
            metadata = route_metadata,
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        throw(ArgumentError("SL_3 Quillen/Murthy route is missing #220 verified global Quillen patch evidence: $(sprint(showerror, err))"))
    end
    verify_quillen_murthy_adapter_consumption(raw_consumption) ||
        throw(ArgumentError("SL_3 Quillen/Murthy route #219/#220 adapter consumption and global patch evidence does not replay"))
    consumption = _polynomial_sl3_quillen_murthy_route_consumption(raw_consumption, route_metadata)
    adapter = _polynomial_quillen_patch_route_adapter(A, consumption.patch)
    replay_metadata = (;
        source = :sl3_quillen_murthy_polynomial_route,
        route_issue_id = "#238",
        route_metadata,
        consumption_replay_metadata = consumption.replay_metadata,
        patch_replay_metadata = adapter.quillen_patch.replay_metadata,
    )
    return (;
        target = adapter.target_matrix,
        route = :quillen_patch,
        context,
        witness_selection = selection,
        local_evidence_provider = provider,
        quillen_consumption = consumption,
        quillen_route_adapter = adapter,
        base_term_policy,
        base_term_factors,
        replay_metadata,
    )
end

function _polynomial_sl3_quillen_murthy_route_core_verification(evidence)
    route_ok = evidence.route == :quillen_patch
    context_ok =
        _verify_sl3_realization_input_context(evidence.context) &&
        evidence.context.matrix == evidence.target &&
        evidence.context.local_form_status == :replayed &&
        evidence.context.evidence_status == :replayable &&
        evidence.context.support_status == :supported &&
        evidence.context.staged_diagnostic.status == :supported
    selection_ok =
        context_ok &&
        _verify_sl3_local_form_witness_selection(evidence.witness_selection) &&
        evidence.witness_selection.context == evidence.context
    provider_ok =
        selection_ok &&
        _verify_sl3_murthy_quillen_local_evidence_provider(evidence.local_evidence_provider) &&
        evidence.local_evidence_provider.context == evidence.context &&
        evidence.local_evidence_provider.witness_selection == evidence.witness_selection &&
        evidence.local_evidence_provider.staged_diagnostic.status == :supported &&
        evidence.local_evidence_provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence
    consumption_ok =
        provider_ok &&
        verify_quillen_murthy_adapter_consumption(evidence.quillen_consumption) &&
        evidence.quillen_consumption.original_input == evidence.target &&
        length(evidence.quillen_consumption.murthy_adapters) == 1 &&
        _same_murthy_quillen_local_adapter_data(
            first(evidence.quillen_consumption.murthy_adapters),
            evidence.local_evidence_provider.murthy_adapter,
        )
    adapter_ok =
        consumption_ok &&
        _verify_polynomial_quillen_patch_route_adapter(evidence.quillen_route_adapter) &&
        evidence.quillen_route_adapter.target_matrix == evidence.target &&
        evidence.quillen_route_adapter.quillen_patch == evidence.quillen_consumption.patch
    base_term_ok =
        adapter_ok &&
        evidence.quillen_route_adapter.quillen_patch.base_term_policy == evidence.base_term_policy &&
        evidence.base_term_policy == :already_handled &&
        isempty(evidence.base_term_factors)
    raw_route_metadata =
        consumption_ok ?
        _polynomial_sl3_quillen_murthy_raw_route_metadata(
            evidence.quillen_consumption.raw_consumption,
        ) :
        nothing
    route_metadata_payload =
        _polynomial_sl3_quillen_murthy_route_metadata_payload(raw_route_metadata)
    expected_route_metadata =
        base_term_ok ?
        _polynomial_sl3_quillen_murthy_route_metadata(
            evidence.context,
            evidence.witness_selection,
            evidence.local_evidence_provider,
            evidence.base_term_policy,
            evidence.base_term_factors,
            route_metadata_payload,
        ) :
        nothing
    route_metadata_ok =
        base_term_ok &&
        raw_route_metadata == expected_route_metadata &&
        hasproperty(evidence.quillen_consumption.replay_metadata, :route_metadata) &&
        evidence.quillen_consumption.replay_metadata.route_metadata == expected_route_metadata &&
        hasproperty(evidence.quillen_consumption.patch.replay_metadata, :metadata) &&
        evidence.quillen_consumption.patch.replay_metadata.metadata == expected_route_metadata &&
        hasproperty(evidence.quillen_route_adapter.quillen_patch.replay_metadata, :metadata) &&
        evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata ==
            expected_route_metadata
    expected_metadata =
        route_metadata_ok ?
        (;
            source = :sl3_quillen_murthy_polynomial_route,
            route_issue_id = "#238",
            route_metadata = expected_route_metadata,
            consumption_replay_metadata = evidence.quillen_consumption.replay_metadata,
            patch_replay_metadata = evidence.quillen_route_adapter.quillen_patch.replay_metadata,
        ) :
        nothing
    replay_metadata_ok = evidence.replay_metadata == expected_metadata
    overall_core_ok =
        route_ok &&
        context_ok &&
        selection_ok &&
        provider_ok &&
        consumption_ok &&
        adapter_ok &&
        base_term_ok &&
        route_metadata_ok &&
        replay_metadata_ok
    return (;
        route_ok,
        context_ok,
        selection_ok,
        provider_ok,
        consumption_ok,
        adapter_ok,
        base_term_ok,
        route_metadata_ok,
        replay_metadata_ok,
        overall_core_ok,
    )
end

function _polynomial_sl3_quillen_murthy_route_verification(evidence)
    core = _polynomial_sl3_quillen_murthy_route_core_verification(evidence)
    stored_verification_ok = evidence.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_polynomial_sl3_quillen_murthy_route_evidence(evidence)::Bool
    try
        return _polynomial_sl3_quillen_murthy_route_verification(evidence).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_sl3_quillen_murthy_route_evidence(A; metadata = (;))
    fields = _polynomial_sl3_quillen_murthy_route_fields(A; metadata)
    raw = PolynomialSL3QuillenMurthyRouteEvidence(values(merge(fields, (; verification = nothing,)))...)
    verification = _polynomial_sl3_quillen_murthy_route_core_verification(raw)
    evidence = PolynomialSL3QuillenMurthyRouteEvidence(values(merge(fields, (; verification,)))...)
    _verify_polynomial_sl3_quillen_murthy_route_evidence(evidence) ||
        error("internal SL_3 Quillen/Murthy route evidence verification failed")
    return evidence
end

function _polynomial_sl3_quillen_murthy_route_certificate(A)
    evidence = _polynomial_sl3_quillen_murthy_route_evidence(A)
    adapter = evidence.quillen_route_adapter
    factors = copy(adapter.global_elementary_factors)
    return _polynomial_route_certificate(
        adapter.target_matrix,
        :quillen_patch,
        factors,
        adapter.product,
        evidence,
        :supported,
    )
end

function _polynomial_sl3_quillen_murthy_route_error(A)
    try
        _polynomial_sl3_quillen_murthy_route_certificate(A)
        return nothing
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return err
    end
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

    if nrows(A) == 3 && length(collect(gens(R))) > 1
        sl3_error = _polynomial_sl3_quillen_murthy_route_error(A)
        if sl3_error === nothing
            return (; error_type = :none, message = "")
        end
        return (;
            error_type = :ArgumentError,
            message = _polynomial_sl3_quillen_murthy_staged_failure_message(sl3_error),
        )
    end

    if nrows(A) > 3
        if allow_recursive_column_peel
            recursive_staged_evidence =
                _polynomial_recursive_column_peel_public_staged_failure_evidence(A)
            if recursive_staged_evidence === nothing
                return (; error_type = :none, message = "")
            elseif recursive_staged_evidence !== _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE
                return recursive_staged_evidence
            end
        end

        staged_error = nothing
        try
            reduce_sln_to_sl3(A)
            return (; error_type = :none, message = "")
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            staged_error = err
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
                !(cert.evidence.final_certificate.evidence isa PolynomialSL3IdentityQuillenRouteEvidence) &&
                _verify_polynomial_column_peel_certificate(cert.evidence) &&
                _polynomial_column_peel_public_mainline_supported(cert.evidence) &&
                _polynomial_route_factor_sequences_equal(cert.factors, cert.evidence.factors)
        elseif cert.route == :quillen_patch
            if cert.evidence isa PolynomialQuillenPatchRouteAdapter
                return cert.evidence.target_matrix == cert.matrix &&
                    _verify_polynomial_quillen_patch_route_adapter(cert.evidence) &&
                    _polynomial_route_factor_sequences_equal(
                        cert.factors,
                        cert.evidence.global_elementary_factors,
                    )
            elseif cert.evidence isa PolynomialSL3QuillenMurthyRouteEvidence
                return cert.evidence.target == cert.matrix &&
                    _verify_polynomial_sl3_quillen_murthy_route_evidence(cert.evidence) &&
                    _polynomial_route_factor_sequences_equal(
                        cert.factors,
                        cert.evidence.quillen_route_adapter.global_elementary_factors,
                    )
            elseif cert.evidence isa PolynomialSL3SuppliedQuillenRouteEvidence
                return cert.evidence.target == cert.matrix &&
                    _verify_polynomial_sl3_supplied_quillen_route_evidence(cert.evidence) &&
                    _polynomial_route_factor_sequences_equal(
                        cert.factors,
                        cert.evidence.quillen_route_adapter.global_elementary_factors,
                    )
            elseif cert.evidence isa PolynomialSL3IdentityQuillenRouteEvidence
                return cert.evidence.target == cert.matrix &&
                    _verify_polynomial_sl3_identity_quillen_route_evidence(cert.evidence) &&
                    _polynomial_route_factor_sequences_equal(
                        cert.factors,
                        cert.evidence.factors,
                    )
            end
            return false
        elseif cert.route == :staged_failure
            hasproperty(cert.evidence, :error_type) &&
                hasproperty(cert.evidence, :message) &&
                cert.evidence.error_type isa Symbol &&
                cert.evidence.message isa AbstractString &&
                !isempty(cert.evidence.message) ||
                return false
            hasproperty(cert.evidence, :allow_recursive_column_peel) ||
                return false
            cert.evidence.allow_recursive_column_peel isa Bool ||
                return false
            allow_recursive_column_peel = cert.evidence.allow_recursive_column_peel
            fresh_evidence = _polynomial_staged_failure_evidence(
                cert.matrix;
                allow_recursive_column_peel,
            )
            fresh_evidence = merge(fresh_evidence, (; allow_recursive_column_peel))
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
