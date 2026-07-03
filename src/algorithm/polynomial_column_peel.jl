struct SLnRecursiveDriverInputContext
    matrix
    base_ring
    coefficient_ring
    dimension::Int
    ring_profile::Symbol
    exact_field_status::Symbol
    determinant
    determinant_status::Symbol
    generators::Tuple
    generator_names::Tuple
    variable_order
    variable_order_status::Symbol
    selected_variable
    selected_variable_index
    selected_variable_status::Symbol
    last_column::Vector
    last_column_profile::NamedTuple
    initial_unimodularity_witness
    initial_unimodularity_witness_status::Symbol
    ecp_witness_metadata
    ecp_evidence_status::Symbol
    final_route_metadata
    final_route_evidence_status::Symbol
    route_provenance_metadata::NamedTuple
    route_provenance_status::Symbol
    catalog_id
    support_classification::Symbol
    staged_reason_code
    staged_diagnostic::NamedTuple
    verification
end

const _SLN_RECURSIVE_DRIVER_STAGED_REASON_CODES = Set([
    :missing_ecp_evidence,
    :missing_final_sl3_route,
    :missing_variable_metadata,
    :unsupported_coefficient_ring,
    :determinant_not_one,
])

function _sln_recursive_driver_extract(data, fields::Tuple)
    data === nothing && return nothing
    for field in fields
        hasproperty(data, field) && return getproperty(data, field)
    end
    return nothing
end

function _sln_recursive_driver_nonempty_identifier(value)::Bool
    value === nothing && return false
    value isa AbstractString && return !isempty(value)
    return true
end

function _sln_recursive_driver_has_replay_payload(metadata)::Bool
    payload = _sln_recursive_driver_extract(
        metadata,
        (
            :replay_steps,
            :replay_metadata,
            :replay_certificate,
            :replay_payload,
            :certificate,
            :verification,
            :patch,
            :route_certificate,
            :factorization_certificate,
        ),
    )
    payload === nothing && return false
    if payload isa Tuple || payload isa AbstractArray || payload isa Set
        return !isempty(payload)
    end
    return true
end

function _sln_recursive_driver_status(metadata)
    status = _sln_recursive_driver_extract(metadata, (:status,))
    return status isa Symbol ? status : nothing
end

function _sln_recursive_driver_normalize_variable(R, variable)
    normalized = hasproperty(variable, :generator) ? getproperty(variable, :generator) : variable
    return _require_substitution_generator(R, normalized)
end

function _sln_recursive_driver_variable_order(R, variable_order)
    variable_order === nothing && return (), :missing
    ring_gens = Tuple(collect(gens(R)))
    variable_order === :auto && return ring_gens, :auto

    normalized = try
        Tuple(_sln_recursive_driver_normalize_variable(R, variable) for variable in variable_order)
    catch err
        err isa InterruptException && rethrow()
        return (), :missing
    end

    isempty(normalized) && return normalized, :missing
    length(unique(normalized)) == length(normalized) || return normalized, :missing
    all(variable -> any(gen -> gen == variable, ring_gens), normalized) || return normalized, :missing
    return normalized, :provided
end

function _sln_recursive_driver_selected_variable(R, selected_variable, variable_order)
    selected_variable === nothing && return nothing, nothing, :missing
    isempty(variable_order) && return nothing, nothing, :missing

    selected = try
        _sln_recursive_driver_normalize_variable(R, selected_variable)
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return nothing, nothing, :missing # COV_EXCL_LINE
    end
    any(variable -> variable == selected, variable_order) || return nothing, nothing, :missing

    selected_index = try
        _ecp_selected_variable_index(R, selected)
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return nothing, nothing, :missing # COV_EXCL_LINE
    end

    return selected, selected_index, :passes
end

function _sln_recursive_driver_determinant_status(A)
    try
        _require_polynomial_sl_determinant(A)
        return :one, det(A)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return :not_one, det(A)
    end
end

function _sln_recursive_driver_last_column_profile(last_column, R)
    entries_in_base_ring = all(entry -> parent(entry) == R, last_column)
    zero_entries = count(entry -> entry == zero(R), last_column)
    last_entry = isempty(last_column) ? nothing : last(last_column)
    target_column = _column_peel_target_column(R, length(last_column))
    return (;
        length = length(last_column),
        entries_in_base_ring,
        zero_entries,
        last_entry_is_one = last_entry == one(R),
        is_target_column = matrix(R, length(last_column), 1, last_column) == target_column,
    )
end

function _sln_recursive_driver_unimodularity_witness(
    last_column,
    R,
    exact_field_status::Symbol,
    determinant_status::Symbol,
)
    determinant_status == :one || return nothing, :missing
    exact_field_status == :supported || return nothing, :unsupported
    return (;
        source = :determinant_one_column_unimodularity,
        ring = R,
        dimension = length(last_column),
        column = Tuple(last_column),
    ), :deduced
end

function _sln_recursive_driver_ecp_status(metadata, column, R, n::Int)
    metadata === nothing && return :missing
    status = _sln_recursive_driver_status(metadata)
    status in (:missing, :absent) && return :missing

    certificate = _sln_recursive_driver_extract(metadata, (:certificate, :ecp_certificate))
    if certificate !== nothing
        certificate isa ECPColumnReductionCertificate || return :missing
        verify_ecp_column_reduction(certificate) || return :missing
        _same_base_ring(certificate.ring, R) || return :missing
        certificate.original_column == column || return :missing
        certificate.final_column == _column_peel_target_column(R, n) || return :missing
        return status == :replayed ? :replayed : :recorded
    end

    replay_id = _sln_recursive_driver_extract(
        metadata,
        (:replay_id, :case_id, :source_case_id, :fixture_id, :mainline_case_id, :upstream_case_id),
    )
    if _sln_recursive_driver_nonempty_identifier(replay_id) ||
            _sln_recursive_driver_has_replay_payload(metadata)
        return :recorded
    end

    return :missing
end

function _sln_recursive_driver_final_route_matrix_ok(matrix_value, R)::Bool
    matrix_value === nothing && return true
    try
        return nrows(matrix_value) == 3 &&
            ncols(matrix_value) == 3 &&
            _same_base_ring(base_ring(matrix_value), R)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sln_recursive_driver_final_route_status(metadata, R)
    metadata === nothing && return :missing
    status = _sln_recursive_driver_status(metadata)
    status in (:missing, :absent) && return :missing

    replay_id = _sln_recursive_driver_extract(
        metadata,
        (:replay_id, :case_id, :final_case_id, :route_id, :fixture_id, :mainline_case_id),
    )
    matrix_value = _sln_recursive_driver_extract(
        metadata,
        (:matrix, :target_matrix, :final_matrix, :final_block),
    )
    has_replay_evidence =
        _sln_recursive_driver_nonempty_identifier(replay_id) ||
        _sln_recursive_driver_has_replay_payload(metadata)

    if status == :replayed &&
            has_replay_evidence &&
            _sln_recursive_driver_final_route_matrix_ok(matrix_value, R)
        return :replayed
    end

    return :recorded
end

function _sln_recursive_driver_route_provenance(route_provenance_metadata, catalog_id)
    if route_provenance_metadata === nothing
        return (; catalog_id = catalog_id), :missing
    end

    provenance = route_provenance_metadata isa NamedTuple ?
        route_provenance_metadata :
        (; value = route_provenance_metadata)
    if provenance == (; catalog_id = catalog_id)
        return provenance, :missing
    end
    return provenance, :recorded
end

function _sln_recursive_driver_staged_diagnostic(
    exact_field_status::Symbol,
    determinant_status::Symbol,
    variable_order_status::Symbol,
    selected_status::Symbol,
    ecp_status::Symbol,
    final_status::Symbol,
    provenance_status::Symbol;
    route_provenance_metadata = (;),
    catalog_id = nothing,
)
    reason_code =
        exact_field_status == :supported ? nothing : :unsupported_coefficient_ring
    reason_code === nothing && determinant_status != :one &&
        (reason_code = :determinant_not_one)
    reason_code === nothing &&
        (variable_order_status == :missing || selected_status != :passes) &&
        (reason_code = :missing_variable_metadata)
    reason_code === nothing && ecp_status != :replayed &&
        (reason_code = :missing_ecp_evidence)
    reason_code === nothing && final_status != :replayed &&
        (reason_code = :missing_final_sl3_route)

    status = reason_code === nothing ? :supported : :staged
    if reason_code !== nothing && !(reason_code in _SLN_RECURSIVE_DRIVER_STAGED_REASON_CODES)
        error("unknown SL_n recursive driver staged reason code") # COV_EXCL_LINE
    end

    message = status == :supported ?
        "SL_n recursive driver input context is replayable" :
        "SL_n recursive driver input context is staged"

    return (;
        status,
        reason_code,
        message,
        exact_field_status,
        determinant_status,
        variable_order_status,
        selected_variable_status = selected_status,
        ecp_evidence_status = ecp_status,
        final_route_evidence_status = final_status,
        route_provenance_status = provenance_status,
        route_provenance_metadata,
        catalog_id,
    )
end

function _sln_recursive_driver_input_context_fields(
    A;
    variable_order = :auto,
    selected_variable = nothing,
    ecp_witness_metadata = nothing,
    final_route_metadata = nothing,
    route_provenance_metadata = nothing,
    catalog_id = nothing,
)
    dimension = _validate_factorization_matrix(A)
    dimension >= 4 ||
        throw(ArgumentError("SL_n recursive driver context requires size at least 4"))
    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    ring_profile == :polynomial ||
        throw(ArgumentError("SL_n recursive driver context requires an ordinary polynomial base ring"))
    coefficient = coefficient_ring(R)
    exact_field_status = _polynomial_exact_field_backed_ring(R) ? :supported : :unsupported
    generators = Tuple(collect(gens(R)))
    generator_names = Tuple(string(generator) for generator in generators)
    normalized_order, variable_order_status =
        _sln_recursive_driver_variable_order(R, variable_order)
    selected, selected_index, selected_status =
        _sln_recursive_driver_selected_variable(R, selected_variable, normalized_order)
    determinant_status, determinant_value = _sln_recursive_driver_determinant_status(A)
    last_column = [A[row, dimension] for row in 1:dimension]
    last_column_profile = _sln_recursive_driver_last_column_profile(last_column, R)
    witness, witness_status =
        _sln_recursive_driver_unimodularity_witness(last_column, R, exact_field_status, determinant_status)
    ecp_status = _sln_recursive_driver_ecp_status(ecp_witness_metadata, last_column, R, dimension)
    final_status = _sln_recursive_driver_final_route_status(final_route_metadata, R)
    provenance, provenance_status =
        _sln_recursive_driver_route_provenance(route_provenance_metadata, catalog_id)
    staged_diagnostic = _sln_recursive_driver_staged_diagnostic(
        exact_field_status,
        determinant_status,
        variable_order_status,
        selected_status,
        ecp_status,
        final_status,
        provenance_status;
        route_provenance_metadata = provenance,
        catalog_id,
    )
    return (; matrix = A, base_ring = R, coefficient_ring = coefficient,
        dimension, ring_profile, exact_field_status, determinant = determinant_value,
        determinant_status, generators, generator_names, variable_order = normalized_order,
        variable_order_status, selected_variable = selected,
        selected_variable_index = selected_index, selected_variable_status = selected_status,
        last_column, last_column_profile, initial_unimodularity_witness = witness,
        initial_unimodularity_witness_status = witness_status,
        ecp_witness_metadata, ecp_evidence_status = ecp_status,
        final_route_metadata, final_route_evidence_status = final_status,
        route_provenance_metadata = provenance, route_provenance_status = provenance_status,
        catalog_id, support_classification = staged_diagnostic.status,
        staged_reason_code = staged_diagnostic.reason_code,
        staged_diagnostic)
end

function _sln_recursive_driver_input_context_core_verification(context)
    recomputed = _sln_recursive_driver_input_context_fields(
        context.matrix;
        variable_order = context.variable_order,
        selected_variable = context.selected_variable,
        ecp_witness_metadata = context.ecp_witness_metadata,
        final_route_metadata = context.final_route_metadata,
        route_provenance_metadata = context.route_provenance_metadata,
        catalog_id = context.catalog_id,
    )

    matrix_ok = context.matrix == recomputed.matrix
    base_ring_ok = context.base_ring == recomputed.base_ring
    coefficient_ring_ok = context.coefficient_ring == recomputed.coefficient_ring
    dimension_ok = context.dimension == recomputed.dimension
    ring_profile_ok = context.ring_profile == recomputed.ring_profile
    exact_field_status_ok = context.exact_field_status == recomputed.exact_field_status
    determinant_ok = context.determinant == recomputed.determinant
    determinant_status_ok = context.determinant_status == recomputed.determinant_status
    generators_ok = context.generators == recomputed.generators
    generator_names_ok = context.generator_names == recomputed.generator_names
    variable_order_ok = context.variable_order == recomputed.variable_order
    variable_order_status_ok = context.variable_order_status == recomputed.variable_order_status
    selected_variable_ok = context.selected_variable == recomputed.selected_variable
    selected_variable_index_ok = context.selected_variable_index == recomputed.selected_variable_index
    selected_variable_status_ok = context.selected_variable_status == recomputed.selected_variable_status
    last_column_ok = context.last_column == recomputed.last_column
    last_column_profile_ok = context.last_column_profile == recomputed.last_column_profile
    initial_unimodularity_witness_ok =
        context.initial_unimodularity_witness == recomputed.initial_unimodularity_witness
    initial_unimodularity_witness_status_ok =
        context.initial_unimodularity_witness_status == recomputed.initial_unimodularity_witness_status
    ecp_witness_metadata_ok = context.ecp_witness_metadata == recomputed.ecp_witness_metadata
    ecp_evidence_status_ok = context.ecp_evidence_status == recomputed.ecp_evidence_status
    final_route_metadata_ok = context.final_route_metadata == recomputed.final_route_metadata
    final_route_evidence_status_ok =
        context.final_route_evidence_status == recomputed.final_route_evidence_status
    route_provenance_metadata_ok =
        context.route_provenance_metadata == recomputed.route_provenance_metadata
    route_provenance_status_ok =
        context.route_provenance_status == recomputed.route_provenance_status
    catalog_id_ok = context.catalog_id == recomputed.catalog_id
    support_classification_ok =
        context.support_classification == recomputed.support_classification
    staged_reason_code_ok = context.staged_reason_code == recomputed.staged_reason_code
    staged_diagnostic_ok = context.staged_diagnostic == recomputed.staged_diagnostic
    overall_core_ok =
        matrix_ok &&
        base_ring_ok &&
        coefficient_ring_ok &&
        dimension_ok &&
        ring_profile_ok &&
        exact_field_status_ok &&
        determinant_ok &&
        determinant_status_ok &&
        generators_ok &&
        generator_names_ok &&
        variable_order_ok &&
        variable_order_status_ok &&
        selected_variable_ok &&
        selected_variable_index_ok &&
        selected_variable_status_ok &&
        last_column_ok &&
        last_column_profile_ok &&
        initial_unimodularity_witness_ok &&
        initial_unimodularity_witness_status_ok &&
        ecp_witness_metadata_ok &&
        ecp_evidence_status_ok &&
        final_route_metadata_ok &&
        final_route_evidence_status_ok &&
        route_provenance_metadata_ok &&
        route_provenance_status_ok &&
        catalog_id_ok &&
        support_classification_ok &&
        staged_reason_code_ok &&
        staged_diagnostic_ok

    return (;
        matrix_ok,
        base_ring_ok,
        coefficient_ring_ok,
        dimension_ok,
        ring_profile_ok,
        exact_field_status_ok,
        determinant_ok,
        determinant_status_ok,
        generators_ok,
        generator_names_ok,
        variable_order_ok,
        variable_order_status_ok,
        selected_variable_ok,
        selected_variable_index_ok,
        selected_variable_status_ok,
        last_column_ok,
        last_column_profile_ok,
        initial_unimodularity_witness_ok,
        initial_unimodularity_witness_status_ok,
        ecp_witness_metadata_ok,
        ecp_evidence_status_ok,
        final_route_metadata_ok,
        final_route_evidence_status_ok,
        route_provenance_metadata_ok,
        route_provenance_status_ok,
        catalog_id_ok,
        support_classification_ok,
        staged_reason_code_ok,
        staged_diagnostic_ok,
        overall_core_ok,
    )
end

function _sln_recursive_driver_input_context_verification(context)
    core = _sln_recursive_driver_input_context_core_verification(context)
    stored_verification_ok = context.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_sln_recursive_driver_input_context(context)::Bool
    try
        return _sln_recursive_driver_input_context_verification(context).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sln_recursive_driver_input_context(A; kwargs...)
    fields = _sln_recursive_driver_input_context_fields(A; kwargs...)
    raw = SLnRecursiveDriverInputContext(values(merge(fields, (; verification = nothing,)))...)
    verification = _sln_recursive_driver_input_context_core_verification(raw)
    checked = SLnRecursiveDriverInputContext(values(merge(fields, (; verification,)))...)
    _verify_sln_recursive_driver_input_context(checked) ||
        error("internal SL_n recursive driver input context verification failed")
    return checked
end

struct PolynomialColumnPeelStep
    dimension::Int
    input_matrix
    last_column::Vector
    left_factors::Vector
    left_certificate
    ecp_evidence
    ecp_route_provenance::NamedTuple
    after_left_matrix
    right_factors::Vector
    right_clearing_coefficients::Tuple
    peeled_matrix
    next_block
    block_embedding_indices::Vector{Int}
    determinant_metadata::NamedTuple
    descent_metadata::NamedTuple
    verification
end

function PolynomialColumnPeelStep(
    dimension::Int,
    input_matrix,
    last_column::Vector,
    left_factors::Vector,
    after_left_matrix,
    right_factors::Vector,
    peeled_matrix,
    next_block,
)
    return _polynomial_column_peel_step_record(
        dimension,
        input_matrix,
        last_column,
        left_factors,
        nothing,
        after_left_matrix,
        right_factors,
        peeled_matrix,
        next_block,
    )
end

function PolynomialColumnPeelStep(
    dimension::Int,
    input_matrix,
    last_column::Vector,
    left_factors::Vector,
    left_certificate,
    after_left_matrix,
    right_factors::Vector,
    peeled_matrix,
    next_block,
)
    return _polynomial_column_peel_step_record(
        dimension,
        input_matrix,
        last_column,
        left_factors,
        left_certificate,
        after_left_matrix,
        right_factors,
        peeled_matrix,
        next_block,
    )
end

struct PolynomialColumnPeelCertificate
    original_matrix
    peel_steps::Vector{PolynomialColumnPeelStep}
    final_block
    final_certificate
    final_factors::Vector
    factors::Vector
    product
    verification
    descent_metadata::NamedTuple
    mainline_support_metadata::NamedTuple
    final_route_provenance::Symbol
end

const _POLYNOMIAL_COLUMN_PEEL_FAST_FINAL_ROUTE_PROVENANCE = :fast_local_sl3
const _POLYNOMIAL_COLUMN_PEEL_BLOCK_FINAL_ROUTE_PROVENANCE = :disjoint_local_blocks
const _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE = :issue184_evidence_backed_sl3
const _POLYNOMIAL_COLUMN_PEEL_UNSUPPORTED_FINAL_ROUTE_PROVENANCE = :unsupported_final_route
const _POLYNOMIAL_COLUMN_PEEL_PUBLIC_ECP_ROUTES = (:embedded_three_block, :unit_entry)

function _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate)::Bool
    try
        final_certificate.route == :quillen_patch || return false
        if final_certificate.evidence isa PolynomialSL3QuillenMurthyRouteEvidence
            _verify_polynomial_sl3_quillen_murthy_route_evidence(final_certificate.evidence) ||
                return false
        elseif final_certificate.evidence isa PolynomialSL3SuppliedQuillenRouteEvidence
            _verify_polynomial_sl3_supplied_quillen_route_evidence(final_certificate.evidence) ||
                return false
        elseif final_certificate.evidence isa PolynomialSL3IdentityQuillenRouteEvidence
            _verify_polynomial_sl3_identity_quillen_route_evidence(final_certificate.evidence) ||
                return false
        else
            return false
        end
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_supported_final_route_ok(final_certificate)::Bool
    final_certificate.route in (:fast_local_sl3, :disjoint_local_blocks) && return true
    return _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate)
end

function _polynomial_column_peel_final_route_provenance(final_certificate)::Symbol
    try
        final_certificate.route == :fast_local_sl3 &&
            return _POLYNOMIAL_COLUMN_PEEL_FAST_FINAL_ROUTE_PROVENANCE
        final_certificate.route == :disjoint_local_blocks &&
            return _POLYNOMIAL_COLUMN_PEEL_BLOCK_FINAL_ROUTE_PROVENANCE
        _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate) &&
            return _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE
    catch err
        err isa InterruptException && rethrow()
    end
    return _POLYNOMIAL_COLUMN_PEEL_UNSUPPORTED_FINAL_ROUTE_PROVENANCE
end

function PolynomialColumnPeelCertificate(
    original_matrix,
    peel_steps::Vector{PolynomialColumnPeelStep},
    final_block,
    final_certificate,
    final_factors::Vector,
    factors::Vector,
    product,
    verification,
)
    final_route_provenance = _polynomial_column_peel_final_route_provenance(final_certificate)
    descent_metadata = _polynomial_column_peel_certificate_descent_metadata(
        peel_steps,
        original_matrix,
        final_block,
        final_route_provenance,
    )
    mainline_support_metadata = _polynomial_column_peel_mainline_support_metadata(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        final_route_provenance,
        descent_metadata,
    )
    return PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
        descent_metadata,
        mainline_support_metadata,
        final_route_provenance,
    )
end

function PolynomialColumnPeelCertificate(
    original_matrix,
    peel_steps::Vector{PolynomialColumnPeelStep},
    final_block,
    final_certificate,
    final_factors::Vector,
    factors::Vector,
    product,
    verification,
    final_route_provenance::Symbol,
)
    descent_metadata = _polynomial_column_peel_certificate_descent_metadata(
        peel_steps,
        original_matrix,
        final_block,
        final_route_provenance,
    )
    mainline_support_metadata = _polynomial_column_peel_mainline_support_metadata(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        final_route_provenance,
        descent_metadata,
    )
    return PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
        descent_metadata,
        mainline_support_metadata,
        final_route_provenance,
    )
end

function _polynomial_column_peel_certificate(A; final_route=nothing)
    _validate_polynomial_column_peel_input(A)
    _validate_polynomial_column_peel_final_route(final_route)

    peel_steps, final_block, final_certificate, final_factors =
        _polynomial_column_peel_recursive(A; final_route=final_route)
    isempty(peel_steps) &&
        throw(ArgumentError("polynomial column-peel certificate requires at least one real peel step"))
    R = base_ring(A)
    factors = _replay_polynomial_column_peel_factors(peel_steps, final_factors, R)
    product = _factor_product(factors, R, nrows(A))

    certificate = PolynomialColumnPeelCertificate(
        A,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        nothing,
    )
    verification = _polynomial_column_peel_core_verification(certificate)
    verification.overall_core_ok || error("internal polynomial column-peel verification failed")
    return PolynomialColumnPeelCertificate(
        certificate.original_matrix,
        certificate.peel_steps,
        certificate.final_block,
        certificate.final_certificate,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        verification,
    )
end

function _verify_polynomial_column_peel_certificate(cert)::Bool
    try
        return _polynomial_column_peel_verification(cert).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_verification(cert)
    core = _polynomial_column_peel_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function _validate_polynomial_column_peel_input(A)
    nrows(A) == ncols(A) ||
        throw(ArgumentError("polynomial column-peel certificate requires a square matrix"))
    nrows(A) >= 3 ||
        throw(ArgumentError("polynomial column-peel certificate requires size at least 3"))
    R = base_ring(A)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("polynomial column-peel certificate does not support Laurent polynomial rings"))
    _factorization_ring_profile(R) == :polynomial ||
        throw(ArgumentError("polynomial column-peel certificate requires an exact ordinary polynomial ring"))
    det(A) == one(R) ||
        throw(ArgumentError("polynomial column-peel certificate requires determinant-one input"))
    return nrows(A)
end

function _validate_polynomial_column_peel_final_route(final_route)
    final_route === nothing && return nothing
    final_route isa Symbol ||
        throw(ArgumentError("polynomial column-peel final route must be a Symbol"))
    final_route in (:fast_local_sl3, :disjoint_local_blocks, :quillen_patch) ||
        throw(ArgumentError("polynomial column-peel final route must be :fast_local_sl3, :disjoint_local_blocks, or :quillen_patch"))
    return nothing
end

function _polynomial_column_peel_recursive(current; final_route=nothing)
    d = nrows(current)
    final_certificate = _polynomial_column_peel_try_final_route(current; final_route=final_route)
    if final_certificate !== nothing
        return PolynomialColumnPeelStep[],
            current,
            final_certificate,
            copy(final_certificate.factors)
    end

    d == 3 &&
        throw(ArgumentError("polynomial column-peel certificate requires a supported final route at size 3"))

    step = _polynomial_column_peel_step(current)
    next_steps, final_block, next_certificate, final_factors =
        _polynomial_column_peel_recursive(step.next_block; final_route=final_route)
    return vcat(PolynomialColumnPeelStep[step], next_steps),
        final_block,
        next_certificate,
        final_factors
end

function _polynomial_column_peel_quillen_final_route_certificate(current)
    for builder in (
            _polynomial_sl3_identity_quillen_route_certificate,
            _polynomial_sl3_supplied_quillen_route_certificate,
        )
        try
            return builder(current)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
        end
    end
    return _polynomial_factorization_route_certificate(
        current;
        allow_recursive_column_peel=false,
    )
end

function _polynomial_column_peel_final_route_matrix_allowed(certificate, current)::Bool
    certificate.matrix != identity_matrix(base_ring(current), nrows(current)) && return true
    return nrows(current) == 3 &&
        certificate.route == :quillen_patch &&
        certificate.evidence isa PolynomialSL3IdentityQuillenRouteEvidence &&
        _verify_polynomial_sl3_identity_quillen_route_evidence(certificate.evidence)
end

function _polynomial_column_peel_try_final_route(current; final_route=nothing)
    candidate_routes =
        final_route === nothing ?
        (:fast_local_sl3, :disjoint_local_blocks, :quillen_patch) :
        (final_route,)

    for route in candidate_routes
        route == :quillen_patch && nrows(current) != 3 && continue
        certificate = try
            if route == :quillen_patch
                _polynomial_column_peel_quillen_final_route_certificate(current)
            else
                _polynomial_factorization_route_certificate(
                    current;
                    route=route,
                    allow_recursive_column_peel=false,
                )
            end
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            nothing
        end

        if certificate !== nothing &&
                certificate.status == :supported &&
                certificate.route == route &&
                _polynomial_column_peel_supported_final_route_ok(certificate) &&
                _polynomial_column_peel_final_route_matrix_allowed(certificate, current)
            return certificate
        end
    end

    return nothing
end

function _polynomial_column_peel_ecp_stage_kind(stage)
    return hasproperty(stage, :kind) ? stage.kind : :unknown
end

function _polynomial_column_peel_ecp_route(evidence)
    evidence isa ECPColumnReductionCertificate || return :unknown
    isempty(evidence.stages) && return :unknown
    terminal = evidence.stages[end]
    if hasproperty(terminal, :route_metadata) &&
            hasproperty(terminal.route_metadata, :route)
        return terminal.route_metadata.route
    end
    return _polynomial_column_peel_ecp_stage_kind(terminal)
end

function _polynomial_column_peel_ecp_route_provenance(evidence)
    if evidence isa ECPColumnReductionCertificate && verify_ecp_column_reduction(evidence)
        return (;
            source = :ecp_column_reduction_certificate,
            verifier = :verify_ecp_column_reduction,
            status = :verified,
            route = _polynomial_column_peel_ecp_route(evidence),
            stage_kinds = tuple((_polynomial_column_peel_ecp_stage_kind(stage) for stage in evidence.stages)...),
            factor_count = length(evidence.factors),
        )
    end
    return (;
        source = :missing_ecp_certificate,
        verifier = :verify_ecp_column_reduction,
        status = :missing,
        route = :unknown,
        stage_kinds = (),
        factor_count = 0,
    )
end

function _polynomial_column_peel_right_clearing_coefficients(after_left, d::Int)
    return tuple((after_left[d, col] for col in 1:(d - 1))...)
end

function _polynomial_column_peel_block_embedding_indices(d::Int)
    return collect(1:(d - 1))
end

function _polynomial_column_peel_determinant_metadata(input_matrix, peeled_matrix, next_block)
    R = base_ring(input_matrix)
    return (;
        input_determinant = det(input_matrix),
        peeled_determinant = det(peeled_matrix),
        next_block_determinant = det(next_block),
        expected_determinant = one(R),
    )
end

function _polynomial_column_peel_descent_metadata(d::Int)
    return (;
        input_dimension = d,
        next_dimension = d - 1,
        dimension_drop = 1,
        route = :polynomial_column_peel,
    )
end

function _polynomial_column_peel_step_record(
    dimension::Int,
    input_matrix,
    last_column::Vector,
    left_factors::Vector,
    left_certificate,
    after_left_matrix,
    right_factors::Vector,
    peeled_matrix,
    next_block;
    require_verified::Bool = false,
)
    ecp_evidence = left_certificate
    ecp_route_provenance = _polynomial_column_peel_ecp_route_provenance(ecp_evidence)
    right_clearing_coefficients =
        _polynomial_column_peel_right_clearing_coefficients(after_left_matrix, dimension)
    block_embedding_indices = _polynomial_column_peel_block_embedding_indices(dimension)
    determinant_metadata =
        _polynomial_column_peel_determinant_metadata(input_matrix, peeled_matrix, next_block)
    descent_metadata = _polynomial_column_peel_descent_metadata(dimension)
    provisional = PolynomialColumnPeelStep(
        dimension,
        input_matrix,
        last_column,
        left_factors,
        left_certificate,
        ecp_evidence,
        ecp_route_provenance,
        after_left_matrix,
        right_factors,
        right_clearing_coefficients,
        peeled_matrix,
        next_block,
        block_embedding_indices,
        determinant_metadata,
        descent_metadata,
        nothing,
    )
    verification = _polynomial_column_peel_step_core_verification(provisional)
    require_verified && !verification.overall_core_ok &&
        throw(ArgumentError("polynomial column-peel step failed exact ECP replay"))
    return PolynomialColumnPeelStep(
        provisional.dimension,
        provisional.input_matrix,
        provisional.last_column,
        provisional.left_factors,
        provisional.left_certificate,
        provisional.ecp_evidence,
        provisional.ecp_route_provenance,
        provisional.after_left_matrix,
        provisional.right_factors,
        provisional.right_clearing_coefficients,
        provisional.peeled_matrix,
        provisional.next_block,
        provisional.block_embedding_indices,
        provisional.determinant_metadata,
        provisional.descent_metadata,
        verification,
    )
end

function _polynomial_column_peel_step(current)
    R = base_ring(current)
    d = nrows(current)
    last_column = [current[row, d] for row in 1:d]
    left_certificate = ecp_column_reduction_certificate(last_column, R)
    left_factors = left_certificate.factors
    left_product = _factor_product(left_factors, R, d)
    recorded_column = matrix(R, d, 1, last_column)
    target_column = _column_peel_target_column(R, d)
    left_product * recorded_column == target_column ||
        throw(ArgumentError("polynomial column-peel left factors failed to send the last column to e_d"))
    after_left = left_product * current
    after_left[:, d:d] == target_column ||
        throw(ArgumentError("polynomial column-peel left product failed to normalize the last column"))
    right_factors = _expected_column_peel_right_factors(after_left, d, R)
    right_product = _factor_product(right_factors, R, d)
    peeled = after_left * right_product
    next_block = matrix(R, [peeled[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
    _is_valid_polynomial_column_peel_step_data(
        d,
        current,
        last_column,
        left_factors,
        after_left,
        right_factors,
        peeled,
        next_block,
    ) || throw(ArgumentError("polynomial column-peel step failed exact replay"))
    return _polynomial_column_peel_step_record(
        d,
        current,
        last_column,
        left_factors,
        left_certificate,
        after_left,
        right_factors,
        peeled,
        next_block,
        require_verified = true,
    )
end

function _polynomial_column_peel_invalid_step_core_verification()
    return (;
        overall_core_ok = false,
        shape_ok = false,
        last_column_ok = false,
        ecp_evidence_ok = false,
        ecp_route_provenance_ok = false,
        left_factors_ok = false,
        after_left_ok = false,
        right_clearing_coefficients_ok = false,
        right_factors_ok = false,
        peeled_matrix_ok = false,
        block_embedding_ok = false,
        next_block_ok = false,
        determinant_metadata_ok = false,
        descent_metadata_ok = false,
    )
end

function _polynomial_column_peel_step_core_verification(step)
    invalid = _polynomial_column_peel_invalid_step_core_verification()
    try
        d = step.dimension
        input_matrix = step.input_matrix
        d >= 4 || return invalid
        nrows(input_matrix) == d && ncols(input_matrix) == d || return invalid
        R = base_ring(input_matrix)
        _is_laurent_polynomial_ring(R) && return invalid
        _factorization_ring_profile(R) == :polynomial || return invalid

        actual_last_column = [input_matrix[row, d] for row in 1:d]
        shape_ok = true
        last_column_ok = step.last_column == actual_last_column
        recorded_column = matrix(R, d, 1, step.last_column)
        target_column = _column_peel_target_column(R, d)

        left_product = _factor_product(step.left_factors, R, d)
        left_factors_ok = left_product * recorded_column == target_column
        after_left_ok =
            left_product * input_matrix == step.after_left_matrix &&
            step.after_left_matrix[:, d:d] == target_column

        evidence = step.ecp_evidence
        ecp_evidence_ok =
            evidence isa ECPColumnReductionCertificate &&
            step.left_certificate == evidence &&
            verify_ecp_column_reduction(evidence) &&
            evidence.original_column == step.last_column &&
            _same_base_ring(evidence.ring, R) &&
            _factor_sequences_equal(evidence.factors, step.left_factors) &&
            evidence.final_column == target_column
        ecp_route_provenance_ok =
            ecp_evidence_ok &&
            step.ecp_route_provenance == _polynomial_column_peel_ecp_route_provenance(evidence)

        expected_coefficients =
            _polynomial_column_peel_right_clearing_coefficients(step.after_left_matrix, d)
        right_clearing_coefficients_ok =
            step.right_clearing_coefficients == expected_coefficients
        expected_right_factors = _expected_column_peel_right_factors(step.after_left_matrix, d, R)
        right_factors_ok = step.right_factors == expected_right_factors
        right_product = _factor_product(step.right_factors, R, d)
        peeled_matrix_ok =
            step.after_left_matrix * right_product == step.peeled_matrix &&
            step.peeled_matrix[d, d] == one(R) &&
            all(step.peeled_matrix[row, d] == zero(R) for row in 1:(d - 1)) &&
            all(step.peeled_matrix[d, col] == zero(R) for col in 1:(d - 1))

        expected_next = matrix(R, [step.peeled_matrix[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
        next_block_ok =
            step.next_block == expected_next &&
            nrows(step.next_block) == d - 1 &&
            ncols(step.next_block) == d - 1 &&
            det(step.next_block) == one(R)
        expected_indices = _polynomial_column_peel_block_embedding_indices(d)
        block_embedding_ok =
            step.block_embedding_indices == expected_indices &&
            step.peeled_matrix == block_embedding(step.next_block, d, step.block_embedding_indices)
        determinant_metadata_ok =
            step.determinant_metadata ==
            _polynomial_column_peel_determinant_metadata(input_matrix, step.peeled_matrix, step.next_block)
        descent_metadata_ok =
            step.descent_metadata == _polynomial_column_peel_descent_metadata(d)

        overall_core_ok =
            shape_ok &&
            last_column_ok &&
            ecp_evidence_ok &&
            ecp_route_provenance_ok &&
            left_factors_ok &&
            after_left_ok &&
            right_clearing_coefficients_ok &&
            right_factors_ok &&
            peeled_matrix_ok &&
            block_embedding_ok &&
            next_block_ok &&
            determinant_metadata_ok &&
            descent_metadata_ok

        return (;
            overall_core_ok,
            shape_ok,
            last_column_ok,
            ecp_evidence_ok,
            ecp_route_provenance_ok,
            left_factors_ok,
            after_left_ok,
            right_clearing_coefficients_ok,
            right_factors_ok,
            peeled_matrix_ok,
            block_embedding_ok,
            next_block_ok,
            determinant_metadata_ok,
            descent_metadata_ok,
        )
    catch err
        err isa InterruptException && rethrow()
        return invalid
    end
end

function _polynomial_column_peel_step_verification(step)
    core = _polynomial_column_peel_step_core_verification(step)
    stored_verification_ok = try
        step.verification == core
    catch err
        err isa InterruptException && rethrow()
        false
    end
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function _is_valid_polynomial_column_peel_step_data(
    d::Int,
    current,
    last_column,
    left_factors,
    after_left,
    right_factors,
    peeled,
    next_block,
)::Bool
    d >= 4 || return false
    nrows(current) == d && ncols(current) == d || return false
    R = base_ring(current)
    _is_laurent_polynomial_ring(R) && return false
    actual_last_column = [current[row, d] for row in 1:d]
    last_column == actual_last_column || return false

    recorded_column = matrix(R, d, 1, last_column)
    target_column = _column_peel_target_column(R, d)
    left_product = try
        _factor_product(left_factors, R, d)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    left_product * current == after_left || return false
    left_product * recorded_column == target_column || return false
    after_left[:, d:d] == target_column || return false

    right_factors == _expected_column_peel_right_factors(after_left, d, R) || return false
    right_product = try
        _factor_product(right_factors, R, d)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    after_left * right_product == peeled || return false
    peeled[d, d] == one(R) || return false
    all(peeled[row, d] == zero(R) for row in 1:(d - 1)) || return false
    all(peeled[d, col] == zero(R) for col in 1:(d - 1)) || return false

    expected_next = matrix(R, [peeled[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
    next_block == expected_next || return false
    return peeled == block_embedding(next_block, d, collect(1:(d - 1)))
end

function _replay_polynomial_column_peel_factors(peel_steps, final_factors, R)
    replayed = collect(final_factors)

    for step in Iterators.reverse(collect(peel_steps))
        replayed = vcat(
            _inverse_elementary_sequence(step.left_factors),
            _embed_upper_left_factors(replayed, R, step.dimension),
            _inverse_elementary_sequence(step.right_factors),
        )
    end

    return replayed
end

function _polynomial_column_peel_certificate_descent_metadata(
    peel_steps,
    original_matrix,
    final_block,
    final_route_provenance::Symbol,
)
    step_dimensions = tuple((step.dimension for step in peel_steps)...)
    next_dimensions = tuple((nrows(step.next_block) for step in peel_steps)...)
    descent_dimensions = tuple(nrows(original_matrix), next_dimensions...)
    expected_step_count = nrows(original_matrix) - nrows(final_block)
    strict_dimension_descent =
        length(peel_steps) == expected_step_count &&
        all(idx -> next_dimensions[idx] == step_dimensions[idx] - 1, eachindex(step_dimensions)) &&
        all(idx -> idx == 1 || step_dimensions[idx] == next_dimensions[idx - 1], eachindex(step_dimensions))
    final_block_is_sl3 = nrows(final_block) == 3 && ncols(final_block) == 3
    return (;
        route=:park_woodburn_recursive_column_peel,
        input_dimension=nrows(original_matrix),
        final_dimension=nrows(final_block),
        step_count=length(peel_steps),
        expected_step_count,
        step_dimensions,
        next_dimensions,
        descent_dimensions,
        strict_dimension_descent,
        final_block_is_sl3,
        final_route_provenance,
    )
end

function _polynomial_column_peel_descent_metadata_ok(cert)::Bool
    try
        hasproperty(cert, :descent_metadata) || return false
        expected = _polynomial_column_peel_certificate_descent_metadata(
            cert.peel_steps,
            cert.original_matrix,
            cert.final_block,
            cert.final_route_provenance,
        )
        return cert.descent_metadata == expected && cert.descent_metadata.strict_dimension_descent
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return false # COV_EXCL_LINE
    end
end

function _polynomial_column_peel_mainline_support_metadata(
    original_matrix,
    peel_steps,
    final_block,
    final_certificate,
    final_factors,
    factors,
    product,
    final_route_provenance::Symbol,
    descent_metadata,
)
    peel_steps_ecp_verified = all(_polynomial_column_peel_public_ecp_peel_step_ok, peel_steps)
    final_route_issue184_ok =
        _polynomial_column_peel_public_final_sl3_evidence_ok(final_certificate)
    factor_replay_ok = _factor_sequences_equal(
        factors,
        _replay_polynomial_column_peel_factors(peel_steps, final_factors, base_ring(original_matrix)),
    )
    product_replay_ok =
        product == _factor_product(factors, base_ring(original_matrix), nrows(original_matrix))
    reconstruction_ok = verify_factorization(original_matrix, factors)
    supported = peel_steps_ecp_verified &&
        descent_metadata.strict_dimension_descent &&
        descent_metadata.final_block_is_sl3 &&
        final_route_issue184_ok &&
        final_route_provenance == _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE &&
        factor_replay_ok &&
        product_replay_ok &&
        reconstruction_ok

    reason_codes = Symbol[]
    peel_steps_ecp_verified || push!(reason_codes, :missing_ecp_peel_evidence)
    descent_metadata.strict_dimension_descent || push!(reason_codes, :non_strict_dimension_descent)
    descent_metadata.final_block_is_sl3 || push!(reason_codes, :final_block_not_sl3)
    (
        final_route_issue184_ok &&
        final_route_provenance == _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE
    ) || push!(reason_codes, :missing_issue184_final_sl3_route)
    factor_replay_ok || push!(reason_codes, :factor_replay_mismatch)
    product_replay_ok || push!(reason_codes, :product_replay_mismatch)
    reconstruction_ok || push!(reason_codes, :factorization_reconstruction_mismatch)

    return (;
        issue_id="#186",
        marker=supported ? :issue186_mainline : :not_issue186_mainline,
        supported,
        reason_codes=tuple(reason_codes...),
        peel_steps_ecp_verified,
        strict_dimension_descent=descent_metadata.strict_dimension_descent,
        final_block_is_sl3=descent_metadata.final_block_is_sl3,
        final_route_issue184_ok,
        final_route_provenance,
        factor_replay_ok,
        product_replay_ok,
        reconstruction_ok,
    )
end

function _polynomial_column_peel_mainline_support_metadata_ok(cert)::Bool
    try
        hasproperty(cert, :descent_metadata) || return false
        hasproperty(cert, :mainline_support_metadata) || return false
        expected = _polynomial_column_peel_mainline_support_metadata(
            cert.original_matrix,
            cert.peel_steps,
            cert.final_block,
            cert.final_certificate,
            cert.final_factors,
            cert.factors,
            cert.product,
            cert.final_route_provenance,
            cert.descent_metadata,
        )
        cert.mainline_support_metadata == expected || return false
        return (cert.mainline_support_metadata.marker == :issue186_mainline) ==
               cert.mainline_support_metadata.supported
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_public_mainline_supported(cert)::Bool
    try
        cert isa PolynomialColumnPeelCertificate || return false
        _verify_polynomial_column_peel_certificate(cert) || return false
        hasproperty(cert, :mainline_support_metadata) || return false
        metadata = cert.mainline_support_metadata
        return hasproperty(metadata, :issue_id) &&
            hasproperty(metadata, :marker) &&
            hasproperty(metadata, :supported) &&
            hasproperty(metadata, :final_route_provenance) &&
            metadata.issue_id == "#186" &&
            metadata.marker == :issue186_mainline &&
            metadata.supported == true &&
            metadata.final_route_provenance ==
                _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return false # COV_EXCL_LINE
    end
end

const _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE =
    :not_polynomial_column_peel_candidate

function _polynomial_column_peel_public_reason_code(cert)::Symbol
    try
        cert isa PolynomialColumnPeelCertificate || return :missing_final_sl3_route
        hasproperty(cert, :mainline_support_metadata) || return :missing_final_sl3_route
        hasproperty(cert.mainline_support_metadata, :reason_codes) || return :missing_final_sl3_route
        reasons = cert.mainline_support_metadata.reason_codes
        :missing_ecp_peel_evidence in reasons && return :missing_ecp_evidence
        :missing_issue184_final_sl3_route in reasons && return :missing_final_sl3_route
        :final_block_not_sl3 in reasons && return :missing_final_sl3_route
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
    end
    return :missing_final_sl3_route
end

function _polynomial_column_peel_public_reason_code(err::ArgumentError)::Symbol
    message = sprint(showerror, err)
    occursin("ECP", message) && return :missing_ecp_evidence
    occursin("left factors", message) && return :missing_ecp_evidence
    occursin("last column", message) && return :missing_ecp_evidence
    occursin("determinant/unit precondition", message) && return :determinant_not_one
    occursin("exact field-backed coefficient ring", message) &&
        return :unsupported_coefficient_ring
    return :missing_final_sl3_route
end

function _polynomial_column_peel_public_reason_code(A, err::ArgumentError)::Symbol
    message = sprint(showerror, err)
    if occursin("supported final route at size 3", message)
        try
            step = _polynomial_column_peel_step(A)
            !_polynomial_column_peel_public_ecp_peel_step_ok(step) &&
                return :missing_ecp_evidence
        catch step_err
            step_err isa InterruptException && rethrow() # COV_EXCL_LINE
        end
    end
    return _polynomial_column_peel_public_reason_code(err)
end

function _polynomial_column_peel_public_staged_message(reason_code::Symbol)
    if reason_code == :missing_ecp_evidence
        return "SL_n recursive column-peel route is staged: missing verified #185/#262 ECP peel evidence"
    elseif reason_code == :missing_final_sl3_route
        return "SL_n recursive column-peel route is staged: missing verified #184/#263 final SL_3 route evidence"
    elseif reason_code == :determinant_not_one
        return "SL_n recursive column-peel route is staged: determinant/unit precondition failed; polynomial inputs must have determinant 1"
    elseif reason_code == :unsupported_coefficient_ring
        return "SL_n recursive column-peel route is staged: recursive public support currently requires an exact field-backed coefficient ring"
    end
    return "SL_n recursive column-peel route is staged: unsupported recursive route evidence"
end

function _polynomial_column_peel_public_staged_failure_evidence(reason_code::Symbol)
    return (;
        error_type = :ArgumentError,
        reason_code,
        message = _polynomial_column_peel_public_staged_message(reason_code),
    )
end

function _polynomial_column_peel_public_non_candidate_error(err::ArgumentError)::Bool
    return occursin(
        "polynomial column-peel certificate requires at least one real peel step",
        sprint(showerror, err),
    )
end

function _polynomial_recursive_column_peel_public_staged_failure_evidence(A)
    try
        evidence = _polynomial_column_peel_certificate(A)
        if evidence.final_certificate.evidence isa PolynomialSL3IdentityQuillenRouteEvidence
            return _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE
        end
        _polynomial_column_peel_public_mainline_supported(evidence) && return nothing
        return _polynomial_column_peel_public_staged_failure_evidence(
            _polynomial_column_peel_public_reason_code(evidence),
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        _polynomial_column_peel_public_non_candidate_error(err) &&
            return _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE
        return _polynomial_column_peel_public_staged_failure_evidence(
            _polynomial_column_peel_public_reason_code(A, err),
        )
    end
end

function _polynomial_column_peel_core_verification(cert)
    preconditions_ok = _polynomial_column_peel_preconditions_ok(cert)
    step_chain_ok = _polynomial_column_peel_step_chain_ok(
        cert.peel_steps,
        cert.original_matrix,
        cert.final_block,
    )
    steps_ok = try
        all(step -> _polynomial_column_peel_step_verification(step).overall_ok, cert.peel_steps)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    left_certificates_ok = try
        all(_polynomial_column_peel_left_certificate_ok, cert.peel_steps)
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        false
    end
    final_certificate_ok = _polynomial_column_peel_final_certificate_ok(cert)
    final_route_provenance_ok = _polynomial_column_peel_final_route_provenance_ok(cert)
    descent_metadata_ok = _polynomial_column_peel_descent_metadata_ok(cert)
    mainline_support_metadata_ok = _polynomial_column_peel_mainline_support_metadata_ok(cert)
    factor_sequence_ok = try
        replayed_factors = _replay_polynomial_column_peel_factors(
            cert.peel_steps,
            cert.final_factors,
            base_ring(cert.original_matrix),
        )
        _factor_sequences_equal(cert.factors, replayed_factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    product_ok = try
        cert.product == _factor_product(cert.factors, base_ring(cert.original_matrix), nrows(cert.original_matrix))
    catch err
        err isa InterruptException && rethrow()
        false
    end
    factors_ok = try
        verify_factorization(cert.original_matrix, cert.factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    overall_core_ok = preconditions_ok && step_chain_ok && steps_ok && left_certificates_ok &&
        final_certificate_ok && final_route_provenance_ok &&
        descent_metadata_ok && mainline_support_metadata_ok &&
        factor_sequence_ok && product_ok && factors_ok
    return (
        overall_core_ok=overall_core_ok,
        preconditions_ok=preconditions_ok,
        step_chain_ok=step_chain_ok,
        steps_ok=steps_ok,
        left_certificates_ok=left_certificates_ok,
        final_certificate_ok=final_certificate_ok,
        final_route_provenance_ok=final_route_provenance_ok,
        descent_metadata_ok=descent_metadata_ok,
        mainline_support_metadata_ok=mainline_support_metadata_ok,
        factor_sequence_ok=factor_sequence_ok,
        product_ok=product_ok,
        factors_ok=factors_ok,
    )
end

function _polynomial_column_peel_final_route_provenance_ok(cert)::Bool
    try
        hasproperty(cert, :final_route_provenance) || return false
        expected = _polynomial_column_peel_final_route_provenance(cert.final_certificate)
        expected == _POLYNOMIAL_COLUMN_PEEL_UNSUPPORTED_FINAL_ROUTE_PROVENANCE && return false
        return cert.final_route_provenance == expected
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_left_certificate_ok(step)::Bool
    try
        replay = _polynomial_column_peel_step_core_verification(step)
        return replay.ecp_evidence_ok && replay.ecp_route_provenance_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_public_ecp_peel_step_ok(step)::Bool
    try
        _polynomial_column_peel_left_certificate_ok(step) || return false
        hasproperty(step, :ecp_route_provenance) || return false
        hasproperty(step.ecp_route_provenance, :route) || return false
        return step.ecp_route_provenance.route in _POLYNOMIAL_COLUMN_PEEL_PUBLIC_ECP_ROUTES
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return false # COV_EXCL_LINE
    end
end

function _polynomial_column_peel_public_final_sl3_evidence_ok(final_certificate)::Bool
    try
        _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate) ||
            return false
        evidence = final_certificate.evidence
        if evidence isa PolynomialSL3SuppliedQuillenRouteEvidence
            supplied = evidence.supplied_evidence
            hasproperty(supplied, :row) || return false
            hasproperty(supplied, :col) || return false
            supplied.row == 1 && supplied.col == 3 && return false
        end
        return true
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        return false # COV_EXCL_LINE
    end
end

function _polynomial_column_peel_preconditions_ok(cert)::Bool
    try
        A = cert.original_matrix
        nrows(A) == ncols(A) || return false
        nrows(A) >= 3 || return false
        isempty(cert.peel_steps) && return false
        R = base_ring(A)
        _is_laurent_polynomial_ring(R) && return false
        _factorization_ring_profile(R) == :polynomial || return false
        det(A) == one(R) || return false
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_step_chain_ok(steps, original_matrix, final_block)::Bool
    collected = collect(steps)
    isempty(collected) && return false
    current = original_matrix

    for step in collected
        step.input_matrix == current || return false
        step.dimension == nrows(current) || return false
        step.dimension == ncols(current) || return false
        current = step.next_block
    end

    current == final_block || return false
    nrows(final_block) == ncols(final_block) || return false
    return nrows(final_block) >= 3
end

function _polynomial_column_peel_final_certificate_ok(cert)::Bool
    try
        final_certificate = cert.final_certificate
        final_certificate isa PolynomialFactorizationRouteCertificate || return false
        final_certificate.status == :supported || return false
        _polynomial_column_peel_supported_final_route_ok(final_certificate) || return false
        final_certificate.matrix == cert.final_block || return false
        _verify_polynomial_factorization_route_certificate(final_certificate) || return false
        _factor_sequences_equal(cert.final_factors, final_certificate.factors) || return false
        cert.final_block == cert.original_matrix && return false
        return verify_factorization(cert.final_block, cert.final_factors)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
