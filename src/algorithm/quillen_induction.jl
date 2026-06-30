function patched_substitution(A, X, r, l::Integer, g)
    l >= 0 || throw(ArgumentError("l must be nonnegative"))

    R = base_ring(A)
    ring_gens = collect(gens(R))
    variable_idx = findfirst(gen -> gen == X, ring_gens)
    variable_idx === nothing && throw(ArgumentError("X must be a generator of the matrix base ring"))

    values = copy(ring_gens)
    values[variable_idx] = ring_gens[variable_idx] + _coerce_into_ring(R, r, "r")^l * _coerce_into_ring(R, g, "g")

    # Oscar's dense matrix constructor consumes entries in column-major order.
    entries = [
        _coerce_into_ring(R, evaluate(A[row, col], values), "patched matrix entry")
        for col in 1:ncols(A), row in 1:nrows(A)
    ]
    return matrix(R, nrows(A), ncols(A), vec(entries))
end

function _quillen_substitute_matrix_scaled_variable(A, X, coefficient)
    R = base_ring(A)
    selected = _require_substitution_generator(R, X)
    scaled_coefficient = _coerce_into_ring(R, coefficient, "substitution coefficient")
    ring_gens = collect(gens(R))
    variable_idx = findfirst(gen -> gen == selected, ring_gens)
    variable_idx === nothing &&
        throw(ArgumentError("X must be a generator of the matrix base ring"))

    values = copy(ring_gens)
    values[variable_idx] = scaled_coefficient * ring_gens[variable_idx]
    entries = [
        _coerce_into_ring(R, evaluate(A[row, col], values), "scaled substitution entry")
        for col in 1:ncols(A), row in 1:nrows(A)
    ]
    return matrix(R, nrows(A), ncols(A), vec(entries))
end

struct QuillenPatchSubstitutionStep
    step_index::Int
    selected_variable
    raw_denominator
    exponent::Int
    powered_denominator
    coverage_multiplier
    sign_convention::Symbol
    previous_coefficient
    next_coefficient
    previous_matrix
    next_matrix
    bracket_target
    replay_metadata
end

struct QuillenPatchSubstitutionChainVerification
    solver_result_ok::Bool
    ring_ok::Bool
    matrix_ok::Bool
    selected_variable_ok::Bool
    sign_convention_ok::Bool
    coefficient_count::Int
    step_count::Int
    bracket_count::Int
    cumulative_coefficients::Vector
    cumulative_coefficients_ok::Bool
    intermediate_matrices::Vector
    intermediate_matrices_ok::Bool
    expected_steps::Vector{QuillenPatchSubstitutionStep}
    steps_ok::Bool
    bracket_matrices::Vector
    bracket_matrices_ok::Bool
    final_coefficient
    final_coefficient_ok::Bool
    base_term
    base_term_ok::Bool
    telescope_product
    telescope_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenDenominatorCoverVerification
    denominator_count::Int
    multiplier_count::Int
    denominators::Vector
    coverage_multipliers::Vector
    parent_ring_ok::Bool
    exact_ring_ok::Bool
    coverage_terms::Vector
    coverage_sum
    coverage_ok::Bool
end

struct QuillenDenominatorCoverCertificate
    ring
    denominators::Vector
    coverage_multipliers::Vector
    coverage_sum
    verification::QuillenDenominatorCoverVerification
end

function _require_quillen_denominator_cover_ring(R)
    (R isa MPolyRing || R isa PolyRing) ||
        throw(ArgumentError("staged Quillen denominator cover certificates require a supported exact ordinary polynomial ring"))

    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("staged Quillen denominator cover certificates require a supported exact ordinary polynomial ring"))
    return R
end

function _normalize_quillen_cover_elements(R, values, label::AbstractString)
    collected = collect(values)
    return [_coerce_into_ring(R, value, label) for value in collected]
end

function _quillen_denominator_cover_verification(R, denominators, coverage_multipliers)
    denominator_count = length(denominators)
    multiplier_count = length(coverage_multipliers)
    denominator_snapshot = copy(denominators)
    multiplier_snapshot = copy(coverage_multipliers)
    exact_ring_ok = Oscar.is_exact_type(typeof(zero(coefficient_ring(R))))
    parent_ring_ok =
        denominator_count == multiplier_count &&
        all(denominator -> parent(denominator) == R, denominators) &&
        all(multiplier -> parent(multiplier) == R, coverage_multipliers)
    coverage_terms = parent_ring_ok ?
        [coverage_multipliers[idx] * denominators[idx] for idx in eachindex(denominators)] :
        Any[]
    coverage_sum = parent_ring_ok ? sum(coverage_terms; init = zero(R)) : zero(R)
    coverage_ok = parent_ring_ok && exact_ring_ok && coverage_sum == one(R)
    return QuillenDenominatorCoverVerification(
        denominator_count,
        multiplier_count,
        denominator_snapshot,
        multiplier_snapshot,
        parent_ring_ok,
        exact_ring_ok,
        coverage_terms,
        coverage_sum,
        coverage_ok,
    )
end

function quillen_denominator_cover_certificate(R, denominators, coverage_multipliers)
    _require_quillen_denominator_cover_ring(R)
    normalized_denominators = _normalize_quillen_cover_elements(R, denominators, "cover denominator")
    normalized_multipliers = _normalize_quillen_cover_elements(R, coverage_multipliers, "cover coverage multiplier")
    isempty(normalized_denominators) &&
        throw(ArgumentError("staged Quillen denominator cover certificates require at least one denominator"))
    length(normalized_denominators) == length(normalized_multipliers) ||
        throw(ArgumentError("staged Quillen denominator cover certificates require matching denominator and multiplier counts"))
    verification = _quillen_denominator_cover_verification(R, normalized_denominators, normalized_multipliers)
    verification.coverage_ok ||
        throw(ArgumentError("staged Quillen denominator cover certificate requires coverage sum to equal one"))
    return QuillenDenominatorCoverCertificate(
        R,
        normalized_denominators,
        normalized_multipliers,
        verification.coverage_sum,
        verification,
    )
end

struct QuillenLocalRealizationCertificate
    original_input
    ring
    size::Int
    selected_variable
    local_certificate::LocalCertificate
    denominator
    coverage_multiplier
    correction::QuillenElementaryCorrection
    factors::Vector
    local_product
    local_correction
    patched_substitution_witness
    witness_metadata
    verification
end

struct QuillenLocalElementaryFactor
    row::Int
    col::Int
    numerator
    denominator
    coverage_multiplier
    provenance
    local_certificate::LocalCertificate
    metadata
end

struct QuillenLocalFactorSequenceVerification
    original_input
    selected_variable
    factor_count::Int
    raw_denominators::Vector
    product_denominator
    normalized_local_contributions::Vector{QuillenLocalContribution}
    normalized_global_elementary_factors::Vector
    local_product
    local_correction
    denominator_data::Vector{QuillenDenominatorData}
    factor_provenance::Vector
    factor_provenance_ok::Bool
    product_denominator_ok::Bool
    normalized_contributions_ok::Bool
    normalized_global_elementary_factors_ok::Bool
    local_product_ok::Bool
    local_correction_ok::Bool
    patched_substitution
    patched_substitution_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenLocalFactorSequenceCertificate
    original_input
    ring
    size::Int
    selected_variable
    factors::Vector{QuillenLocalElementaryFactor}
    raw_denominators::Vector
    product_denominator
    local_product
    local_correction
    normalized_local_contributions::Vector{QuillenLocalContribution}
    normalized_global_elementary_factors::Vector
    patched_substitution_witness
    chain_witness
    witness_metadata
    replay_metadata
    verification::QuillenLocalFactorSequenceVerification
end

struct QuillenLocalDenominatorSupport
    local_index::Int
    support_denominator
    support_kind::Symbol
    factor_denominators::Vector
    factor_entries::Vector{QuillenLocalElementaryFactor}
    factor_provenance::Vector
    replayed_denominator
    replay_equality
    replay_ok::Bool
end

struct MurthyQuillenLocalAdapter
    original_input
    ring
    size::Int
    selected_variable
    murthy_certificate::SL3LocalRealizationCertificate
    local_factor_replay::SL3LocalElementaryFactorReplay
    mode::Symbol
    materialized_factors
    local_product
    local_correction
    quillen_factor_sequence
    quillen_local_certificate
    witness_metadata
    replay_metadata
    verification
end

struct QuillenDenominatorCoverCandidateVerification
    local_count::Int
    raw_denominators::Vector
    local_certificates_ok::Bool
    same_original_input_ok::Bool
    same_ring_ok::Bool
    same_size_ok::Bool
    same_selected_variable_ok::Bool
    local_supports::Vector{QuillenLocalDenominatorSupport}
    local_supports_ok::Bool
    raw_denominators_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenDenominatorCoverCandidate
    original_input
    ring
    size::Int
    selected_variable
    local_certificates::Vector{QuillenLocalFactorSequenceCertificate}
    raw_denominators::Vector
    local_supports::Vector{QuillenLocalDenominatorSupport}
    replay_metadata
    verification::QuillenDenominatorCoverCandidateVerification
end

struct QuillenDenominatorCoverSolverVerification
    raw_denominator_count::Int
    multiplier_count::Int
    raw_denominators::Vector
    exponent::Int
    powered_denominators::Vector
    coverage_multipliers::Vector
    parent_ring_ok::Bool
    exact_ring_ok::Bool
    exponent_ok::Bool
    source_candidate_ok::Bool
    coverage_terms::Vector
    coverage_sum
    coverage_ok::Bool
    cover_certificate_ok::Bool
    cover_certificate_matches::Bool
    overall_ok::Bool
end

struct QuillenDenominatorCoverSolverResult
    source_candidate
    ring
    raw_denominators::Vector
    exponent::Int
    powered_denominators::Vector
    coverage_multipliers::Vector
    coverage_terms::Vector
    coverage_sum
    cover_certificate::QuillenDenominatorCoverCertificate
    verification::QuillenDenominatorCoverSolverVerification
end

struct QuillenPatchSubstitutionChain
    original_matrix
    ring
    size::Int
    selected_variable
    sign_convention::Symbol
    solver_result::QuillenDenominatorCoverSolverResult
    cumulative_coefficients::Vector
    intermediate_matrices::Vector
    steps::Vector{QuillenPatchSubstitutionStep}
    bracket_matrices::Vector
    base_term
    metadata
    replay_metadata
    verification::QuillenPatchSubstitutionChainVerification
end

struct QuillenLocalContributionNormalizationVerification
    local_certificate_ok::Bool
    cover_certificate_ok::Bool
    original_input_ok::Bool
    selected_variable_ok::Bool
    cover_ring_ok::Bool
    cover_index_ok::Bool
    cover_pair
    cover_pair_ok::Bool
    patched_substitution
    patched_substitution_ok::Bool
    local_product
    local_product_ok::Bool
    local_correction
    local_correction_ok::Bool
    local_contribution_ok::Bool
    weighted_global_elementary_factors::Vector
    weighted_global_elementary_factors_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenLocalContributionNormalization
    local_certificate::QuillenLocalRealizationCertificate
    cover_certificate::QuillenDenominatorCoverCertificate
    original_input
    selected_variable
    denominator
    coverage_multiplier
    cover_index::Int
    patched_substitution_witness
    patched_substitution
    local_product
    local_correction
    local_contribution::QuillenLocalContribution
    weighted_global_elementary_factors::Vector
    replay_metadata
    verification::QuillenLocalContributionNormalizationVerification
end

function _same_quillen_local_contribution(
    left::QuillenLocalContribution,
    right::QuillenLocalContribution,
)::Bool
    return left.certificate.indices == right.certificate.indices &&
           left.certificate.denominators == right.certificate.denominators &&
           left.denominator == right.denominator &&
           left.coverage_multiplier == right.coverage_multiplier &&
           left.correction == right.correction
end

function _quillen_cover_pair_index(
    cover::QuillenDenominatorCoverCertificate,
    denominator,
    coverage_multiplier,
)
    matches = Int[]
    for idx in eachindex(cover.denominators)
        if cover.denominators[idx] == denominator &&
           cover.coverage_multipliers[idx] == coverage_multiplier
            push!(matches, idx)
        end
    end
    isempty(matches) &&
        throw(ArgumentError("local contribution denominator and coverage multiplier must match an exact cover pair"))
    length(matches) == 1 ||
        throw(ArgumentError("local contribution denominator and coverage multiplier match multiple cover pairs"))
    return only(matches)
end

function _quillen_normalization_replay_metadata(certificate, cover, cover_index::Int)
    return (;
        local_witness_metadata = certificate.witness_metadata,
        cover_index = cover_index,
        cover_denominator_count = length(cover.denominators),
        cover_coverage_sum = cover.coverage_sum,
    )
end

function _same_quillen_normalization_verification(
    left::QuillenLocalContributionNormalizationVerification,
    right::QuillenLocalContributionNormalizationVerification,
)::Bool
    return left.local_certificate_ok == right.local_certificate_ok &&
           left.cover_certificate_ok == right.cover_certificate_ok &&
           left.original_input_ok == right.original_input_ok &&
           left.selected_variable_ok == right.selected_variable_ok &&
           left.cover_ring_ok == right.cover_ring_ok &&
           left.cover_index_ok == right.cover_index_ok &&
           left.cover_pair == right.cover_pair &&
           left.cover_pair_ok == right.cover_pair_ok &&
           left.patched_substitution == right.patched_substitution &&
           left.patched_substitution_ok == right.patched_substitution_ok &&
           left.local_product == right.local_product &&
           left.local_product_ok == right.local_product_ok &&
           left.local_correction == right.local_correction &&
           left.local_correction_ok == right.local_correction_ok &&
           left.local_contribution_ok == right.local_contribution_ok &&
           _same_quillen_factors(
               left.weighted_global_elementary_factors,
               right.weighted_global_elementary_factors,
           ) &&
           left.weighted_global_elementary_factors_ok ==
               right.weighted_global_elementary_factors_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _quillen_local_input_ring_size(original_input; ring = nothing, size = nothing)
    if original_input isa QuillenElementaryCorrection
        ring === nothing &&
            throw(ArgumentError("ring is required when original input is an elementary correction"))
        size === nothing &&
            throw(ArgumentError("size is required when original input is an elementary correction"))
        R = _require_supported_quillen_ring(ring)
        size isa Integer && size >= 2 ||
            throw(ArgumentError("certificate size must be at least 2"))
        return R, Int(size)
    end

    matrix_size = _require_square_matrix(original_input, "original input")
    R = _require_supported_quillen_ring(base_ring(original_input))
    return R, matrix_size
end

function _quillen_local_require_factor_matrix(factor, R, n::Int, label::AbstractString)
    nrows(factor) == n && ncols(factor) == n ||
        throw(ArgumentError("$(label) must be a square matrix of certificate size"))
    base_ring(factor) == R ||
        throw(ArgumentError("$(label) must be defined over the certificate ring"))
    return factor
end

function _quillen_local_factor_vector(factors, local_correction, R, n::Int)
    if factors === nothing
        local_correction === nothing &&
            throw(ArgumentError("either factors or local_correction must be supplied"))
        return [_quillen_local_require_factor_matrix(local_correction, R, n, "local correction")]
    end
    collected = collect(factors)
    isempty(collected) && throw(ArgumentError("local factors must be nonempty"))
    return [
        _quillen_local_require_factor_matrix(factor, R, n, "local factor")
        for factor in collected
    ]
end

function _quillen_local_patched_witness_summary(witness, R, n::Int, selected_variable)
    witness === nothing && return (; present = false, ok = true)
    for field in (:matrix, :variable, :denominator, :exponent, :shift, :expected_matrix)
        hasproperty(witness, field) ||
            throw(ArgumentError("patched-substitution witness missing field $(field)"))
    end
    _quillen_local_require_factor_matrix(
        witness.matrix,
        R,
        n,
        "patched-substitution witness matrix",
    )
    _quillen_local_require_factor_matrix(
        witness.expected_matrix,
        R,
        n,
        "patched-substitution witness expected matrix",
    )
    witness_variable = _require_substitution_generator(R, witness.variable)
    denominator = _coerce_into_ring(R, witness.denominator, "patched-substitution witness denominator")
    shift = _coerce_into_ring(R, witness.shift, "patched-substitution witness shift")
    actual = patched_substitution(
        witness.matrix,
        witness_variable,
        denominator,
        witness.exponent,
        shift,
    )
    ok = witness_variable == selected_variable && actual == witness.expected_matrix
    return (;
        present = true,
        variable = witness_variable,
        denominator = denominator,
        exponent = witness.exponent,
        shift = shift,
        actual_matrix = actual,
        expected_matrix = witness.expected_matrix,
        ok = ok,
    )
end

function _quillen_local_certificate_replay_summary(
    certificate::QuillenLocalRealizationCertificate,
)
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    n >= 2 || throw(ArgumentError("certificate size must be at least 2"))
    selected_variable = _require_substitution_generator(R, certificate.selected_variable)

    original_input = if certificate.original_input isa QuillenElementaryCorrection
        certificate.original_input == certificate.correction ||
            throw(ArgumentError("original elementary correction must match recorded correction"))
        certificate.original_input
    else
        _quillen_local_require_factor_matrix(certificate.original_input, R, n, "original input")
    end

    normalized = _normalize_quillen_contribution(
        QuillenLocalContribution(
            certificate.local_certificate,
            certificate.denominator,
            certificate.coverage_multiplier,
            certificate.correction,
        ),
        R,
        n,
    )
    expected_factor = _quillen_factors(R, n, [normalized])[1]
    factors = [
        _quillen_local_require_factor_matrix(factor, R, n, "local factor")
        for factor in certificate.factors
    ]
    local_product = _quillen_product(R, n, factors)
    local_correction = _quillen_local_require_factor_matrix(
        certificate.local_correction,
        R,
        n,
        "local correction",
    )
    witness = _quillen_local_patched_witness_summary(
        certificate.patched_substitution_witness,
        R,
        n,
        selected_variable,
    )
    denominator_ok = normalized.denominator == certificate.denominator
    correction_ok = expected_factor == local_correction
    factors_ok = local_product == local_correction
    stored_product_ok = certificate.local_product == local_product
    overall_ok =
        denominator_ok && correction_ok && factors_ok && stored_product_ok && witness.ok
    return (;
        original_input = original_input,
        selected_variable = selected_variable,
        denominator = normalized.denominator,
        coverage_multiplier = normalized.coverage_multiplier,
        expected_local_correction = expected_factor,
        local_product = local_product,
        local_correction = local_correction,
        denominator_ok = denominator_ok,
        local_correction_ok = correction_ok,
        factors_ok = factors_ok,
        stored_product_ok = stored_product_ok,
        patched_substitution = witness,
        patched_substitution_ok = witness.ok,
        witness_metadata = certificate.witness_metadata,
        overall_ok = overall_ok,
    )
end

function quillen_local_realization_certificate(
    original_input,
    selected_variable;
    local_certificate::LocalCertificate,
    denominator,
    coverage_multiplier,
    correction::QuillenElementaryCorrection,
    factors = nothing,
    local_correction = nothing,
    patched_substitution_witness = nothing,
    witness_metadata = (;),
    ring = nothing,
    size = nothing,
)
    R, n = _quillen_local_input_ring_size(original_input; ring, size)
    selected = _require_substitution_generator(R, selected_variable)
    normalized = _normalize_quillen_contribution(
        QuillenLocalContribution(
            local_certificate,
            denominator,
            coverage_multiplier,
            correction,
        ),
        R,
        n,
    )
    normalized_factors = _quillen_local_factor_vector(factors, local_correction, R, n)
    product = _quillen_product(R, n, normalized_factors)
    recorded_correction = local_correction === nothing ?
        product :
        _quillen_local_require_factor_matrix(local_correction, R, n, "local correction")
    provisional = QuillenLocalRealizationCertificate(
        original_input,
        R,
        n,
        selected,
        normalized.certificate,
        normalized.denominator,
        normalized.coverage_multiplier,
        normalized.correction,
        normalized_factors,
        product,
        recorded_correction,
        patched_substitution_witness,
        witness_metadata,
        nothing,
    )
    verification = _quillen_local_certificate_replay_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen local realization certificate data does not replay"))
    return QuillenLocalRealizationCertificate(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.local_certificate,
        provisional.denominator,
        provisional.coverage_multiplier,
        provisional.correction,
        provisional.factors,
        provisional.local_product,
        provisional.local_correction,
        provisional.patched_substitution_witness,
        provisional.witness_metadata,
        verification,
    )
end

function replay_quillen_denominator_cover(certificate::QuillenDenominatorCoverCertificate)
    _require_quillen_denominator_cover_ring(certificate.ring)
    return _quillen_denominator_cover_verification(
        certificate.ring,
        certificate.denominators,
        certificate.coverage_multipliers,
    )
end

function _same_quillen_denominator_cover_verification(
    left::QuillenDenominatorCoverVerification,
    right::QuillenDenominatorCoverVerification,
)::Bool
    return left.denominator_count == right.denominator_count &&
           left.multiplier_count == right.multiplier_count &&
           left.denominators == right.denominators &&
           left.coverage_multipliers == right.coverage_multipliers &&
           left.parent_ring_ok == right.parent_ring_ok &&
           left.exact_ring_ok == right.exact_ring_ok &&
           left.coverage_terms == right.coverage_terms &&
           left.coverage_sum == right.coverage_sum &&
           left.coverage_ok == right.coverage_ok
end

function verify_quillen_denominator_cover(certificate::QuillenDenominatorCoverCertificate)::Bool
    try
        replay = replay_quillen_denominator_cover(certificate)
        return certificate.coverage_sum == replay.coverage_sum &&
               _same_quillen_denominator_cover_verification(certificate.verification, replay) &&
               replay.parent_ring_ok &&
               replay.exact_ring_ok &&
               replay.coverage_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function verify_quillen_local_certificate(certificate)::Bool
    try
        replay = _quillen_local_certificate_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function QuillenLocalElementaryFactor(
    row::Int,
    col::Int,
    numerator,
    denominator,
    coverage_multiplier,
    local_certificate::LocalCertificate,
    provenance,
    metadata,
)
    return QuillenLocalElementaryFactor(
        row,
        col,
        numerator,
        denominator,
        coverage_multiplier,
        provenance,
        local_certificate,
        metadata,
    )
end

function _quillen_local_sequence_factor_field(factor, field::Symbol)
    hasproperty(factor, field) ||
        throw(ArgumentError("local elementary factor missing field $(field)"))
    return getproperty(factor, field)
end

function _quillen_local_sequence_require_nonempty(value, label::AbstractString)
    value === nothing && throw(ArgumentError("$(label) must be nonempty"))
    if value isa NamedTuple
        isempty(keys(value)) && throw(ArgumentError("$(label) must be nonempty"))
        return value
    end
    applicable(isempty, value) && isempty(value) &&
        throw(ArgumentError("$(label) must be nonempty"))
    return value
end

function _quillen_local_sequence_original_input(original_input, R, n::Int)
    if original_input isa QuillenElementaryCorrection
        row, col = _require_elementary_indices(n, original_input.row, original_input.col)
        entry = _coerce_into_ring(R, original_input.entry, "original input correction entry")
        return QuillenElementaryCorrection(row, col, entry)
    end
    return _quillen_local_require_factor_matrix(original_input, R, n, "original input")
end

function _quillen_local_sequence_factor(raw_factor, R, n::Int, index::Int)
    row = Int(_quillen_local_sequence_factor_field(raw_factor, :row))
    col = Int(_quillen_local_sequence_factor_field(raw_factor, :col))
    row, col = _require_elementary_indices(n, row, col)
    numerator = _coerce_into_ring(
        R,
        _quillen_local_sequence_factor_field(raw_factor, :numerator),
        "local elementary factor numerator",
    )
    denominator = _coerce_into_ring(
        R,
        _quillen_local_sequence_factor_field(raw_factor, :denominator),
        "local elementary factor denominator",
    )
    coverage_multiplier = _coerce_into_ring(
        R,
        _quillen_local_sequence_factor_field(raw_factor, :coverage_multiplier),
        "local elementary factor coverage multiplier",
    )
    provenance = _quillen_local_sequence_require_nonempty(
        _quillen_local_sequence_factor_field(raw_factor, :provenance),
        "local elementary factor provenance",
    )
    local_certificate = hasproperty(raw_factor, :local_certificate) &&
                        getproperty(raw_factor, :local_certificate) !== nothing ?
        getproperty(raw_factor, :local_certificate) :
        LocalCertificate([row, col], [denominator, denominator])
    local_certificate isa LocalCertificate ||
        throw(ArgumentError("local elementary factor local certificate must be a LocalCertificate"))
    metadata = hasproperty(raw_factor, :metadata) ? getproperty(raw_factor, :metadata) : (; factor_index = index)
    return QuillenLocalElementaryFactor(
        row,
        col,
        numerator,
        denominator,
        coverage_multiplier,
        provenance,
        local_certificate,
        metadata,
    )
end

function _quillen_local_sequence_contributions(factors, R, n::Int)
    return [
        _normalize_quillen_contribution(
            QuillenLocalContribution(
                factor.local_certificate,
                factor.denominator,
                factor.coverage_multiplier,
                QuillenElementaryCorrection(factor.row, factor.col, factor.numerator),
            ),
            R,
            n,
        )
        for factor in factors
    ]
end

function _quillen_local_sequence_provenance_field(provenance, field::Symbol)
    if hasproperty(provenance, field)
        return true, getproperty(provenance, field)
    end
    if applicable(haskey, provenance, field) && haskey(provenance, field)
        return true, provenance[field]
    end
    return false, nothing
end

function _quillen_local_sequence_provenance_matches_position(
    provenance,
    index::Int,
    field::Symbol,
)::Bool
    present, value = _quillen_local_sequence_provenance_field(provenance, field)
    present || return true
    value isa Integer || return false
    return Int(value) == index
end

function _quillen_local_sequence_factor_provenance_ok(provenance, index::Int)::Bool
    _quillen_local_sequence_require_nonempty(provenance, "local elementary factor provenance")
    return _quillen_local_sequence_provenance_matches_position(
        provenance,
        index,
        :factor_index,
    ) &&
           _quillen_local_sequence_provenance_matches_position(
               provenance,
               index,
               :sequence_index,
           ) &&
           _quillen_local_sequence_provenance_matches_position(
               provenance,
               index,
               :local_index,
           )
end

function _same_quillen_local_contributions(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _same_quillen_local_contribution(left[idx], right[idx]) || return false
    end
    return true
end

function _quillen_local_sequence_factor_provenance(factors)
    return [factor.provenance for factor in factors]
end

function _same_quillen_local_elementary_factor(
    left::QuillenLocalElementaryFactor,
    right::QuillenLocalElementaryFactor,
)::Bool
    return left.row == right.row &&
           left.col == right.col &&
           left.numerator == right.numerator &&
           left.denominator == right.denominator &&
           left.coverage_multiplier == right.coverage_multiplier &&
           left.provenance == right.provenance &&
           left.local_certificate.indices == right.local_certificate.indices &&
           left.local_certificate.denominators == right.local_certificate.denominators &&
           left.metadata == right.metadata
end

function _same_quillen_local_elementary_factors(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _same_quillen_local_elementary_factor(left[idx], right[idx]) || return false
    end
    return true
end

function _quillen_local_denominator_support(
    certificate::QuillenLocalFactorSequenceCertificate,
    local_index::Int,
)
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    factors = [
        _quillen_local_sequence_factor(raw_factor, R, n, index)
        for (index, raw_factor) in enumerate(certificate.factors)
    ]
    factor_denominators = [factor.denominator for factor in factors]
    factor_provenance = _quillen_local_sequence_factor_provenance(factors)
    replayed_denominator = prod(
        (factor.denominator for factor in factors);
        init = one(certificate.ring),
    )
    replay_equality =
        certificate.raw_denominators == factor_denominators &&
        certificate.product_denominator == replayed_denominator
    return QuillenLocalDenominatorSupport(
        Int(local_index),
        certificate.product_denominator,
        :product,
        factor_denominators,
        factors,
        factor_provenance,
        replayed_denominator,
        replay_equality,
        replay_equality,
    )
end

function _same_quillen_local_denominator_support(
    left::QuillenLocalDenominatorSupport,
    right::QuillenLocalDenominatorSupport,
)::Bool
    return left.local_index == right.local_index &&
           left.support_denominator == right.support_denominator &&
           left.support_kind == right.support_kind &&
           left.factor_denominators == right.factor_denominators &&
           _same_quillen_local_elementary_factors(left.factor_entries, right.factor_entries) &&
           left.factor_provenance == right.factor_provenance &&
           left.replayed_denominator == right.replayed_denominator &&
           left.replay_equality == right.replay_equality &&
           left.replay_ok == right.replay_ok
end

function _same_quillen_local_denominator_supports(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _same_quillen_local_denominator_support(left[idx], right[idx]) || return false
    end
    return true
end

function _quillen_denominator_cover_candidate_replay_metadata(
    local_certificates::Vector{QuillenLocalFactorSequenceCertificate},
    local_supports::Vector{QuillenLocalDenominatorSupport},
    raw_denominators,
)
    return (;
        local_count = length(local_certificates),
        raw_denominators = collect(raw_denominators),
        support_denominators = [support.support_denominator for support in local_supports],
        support_kinds = [support.support_kind for support in local_supports],
        local_sequence_replay_metadata = [certificate.replay_metadata for certificate in local_certificates],
    )
end

function _same_quillen_denominator_cover_candidate_verification(
    left::QuillenDenominatorCoverCandidateVerification,
    right::QuillenDenominatorCoverCandidateVerification,
)::Bool
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
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _quillen_local_sequence_replay_metadata(
    certificate::QuillenLocalFactorSequenceCertificate,
    denominator_data,
)
    return (;
        factor_count = length(certificate.factors),
        raw_denominators = [factor.denominator for factor in certificate.factors],
        denominator_data = denominator_data,
        factor_provenance = _quillen_local_sequence_factor_provenance(certificate.factors),
        factor_metadata = [factor.metadata for factor in certificate.factors],
        factor_local_certificates = [factor.local_certificate for factor in certificate.factors],
        has_patched_substitution_witness = certificate.patched_substitution_witness !== nothing,
        has_chain_witness = certificate.chain_witness !== nothing,
        chain_witness = certificate.chain_witness,
        witness_metadata = certificate.witness_metadata,
    )
end

function _same_quillen_local_factor_sequence_verification(
    left::QuillenLocalFactorSequenceVerification,
    right::QuillenLocalFactorSequenceVerification,
)::Bool
    return left.original_input == right.original_input &&
           left.selected_variable == right.selected_variable &&
           left.factor_count == right.factor_count &&
           left.raw_denominators == right.raw_denominators &&
           left.product_denominator == right.product_denominator &&
           _same_quillen_local_contributions(
               left.normalized_local_contributions,
               right.normalized_local_contributions,
           ) &&
           _same_quillen_factors(
               left.normalized_global_elementary_factors,
               right.normalized_global_elementary_factors,
           ) &&
           left.local_product == right.local_product &&
           left.local_correction == right.local_correction &&
           _same_quillen_denominator_data(left.denominator_data, right.denominator_data) &&
           left.factor_provenance == right.factor_provenance &&
           left.factor_provenance_ok == right.factor_provenance_ok &&
           left.product_denominator_ok == right.product_denominator_ok &&
           left.normalized_contributions_ok == right.normalized_contributions_ok &&
           left.normalized_global_elementary_factors_ok ==
               right.normalized_global_elementary_factors_ok &&
           left.local_product_ok == right.local_product_ok &&
           left.local_correction_ok == right.local_correction_ok &&
           left.patched_substitution == right.patched_substitution &&
           left.patched_substitution_ok == right.patched_substitution_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function replay_quillen_local_factor_sequence(
    certificate::QuillenLocalFactorSequenceCertificate,
)
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    n >= 2 || throw(ArgumentError("certificate size must be at least 2"))
    selected_variable = _require_substitution_generator(R, certificate.selected_variable)
    original_input = _quillen_local_sequence_original_input(certificate.original_input, R, n)
    factors = [
        _quillen_local_sequence_factor(raw_factor, R, n, index)
        for (index, raw_factor) in enumerate(certificate.factors)
    ]
    raw_denominators = [factor.denominator for factor in factors]
    product_denominator = prod(raw_denominators; init = one(R))
    normalized_local_contributions = _quillen_local_sequence_contributions(factors, R, n)
    normalized_global_elementary_factors = _quillen_factors(R, n, normalized_local_contributions)
    local_product = _quillen_product(R, n, normalized_global_elementary_factors)
    local_correction = _quillen_local_require_factor_matrix(
        certificate.local_correction,
        R,
        n,
        "local correction",
    )
    denominator_data = _quillen_denominator_data(normalized_local_contributions)
    factor_provenance = _quillen_local_sequence_factor_provenance(factors)
    factor_provenance_ok = all(enumerate(factor_provenance)) do (index, provenance)
        try
            _quillen_local_sequence_factor_provenance_ok(provenance, index)
        catch err
            err isa InterruptException && rethrow()
            false
        end
    end
    product_denominator_ok =
        certificate.raw_denominators == raw_denominators &&
        certificate.product_denominator == product_denominator
    normalized_contributions_ok = _same_quillen_local_contributions(
        certificate.normalized_local_contributions,
        normalized_local_contributions,
    )
    normalized_global_elementary_factors_ok = _same_quillen_factors(
        certificate.normalized_global_elementary_factors,
        normalized_global_elementary_factors,
    )
    local_product_ok = certificate.local_product == local_product
    local_correction_ok = local_correction == local_product
    patched_substitution = _quillen_local_patched_witness_summary(
        certificate.patched_substitution_witness,
        R,
        n,
        selected_variable,
    )
    patched_substitution_ok = patched_substitution.ok
    replay_metadata = _quillen_local_sequence_replay_metadata(certificate, denominator_data)
    replay_metadata_ok = certificate.replay_metadata == replay_metadata
    overall_ok =
        factor_provenance_ok &&
        product_denominator_ok &&
        normalized_contributions_ok &&
        normalized_global_elementary_factors_ok &&
        local_product_ok &&
        local_correction_ok &&
        patched_substitution_ok &&
        replay_metadata_ok

    return QuillenLocalFactorSequenceVerification(
        original_input,
        selected_variable,
        length(factors),
        raw_denominators,
        product_denominator,
        normalized_local_contributions,
        normalized_global_elementary_factors,
        local_product,
        local_correction,
        denominator_data,
        factor_provenance,
        factor_provenance_ok,
        product_denominator_ok,
        normalized_contributions_ok,
        normalized_global_elementary_factors_ok,
        local_product_ok,
        local_correction_ok,
        patched_substitution,
        patched_substitution_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function QuillenLocalFactorSequenceCertificate(
    ring,
    size::Int,
    selected_variable,
    factors,
    raw_denominators,
    product_denominator,
    normalized_global_elementary_factors,
    local_product,
    local_correction,
    verification,
)
    R = _require_supported_quillen_ring(ring)
    n = Int(size)
    normalized_factors = [
        _quillen_local_sequence_factor(raw_factor, R, n, index)
        for (index, raw_factor) in enumerate(collect(factors))
    ]
    normalized_local_contributions = try
        _quillen_local_sequence_contributions(normalized_factors, R, n)
    catch err
        err isa InterruptException && rethrow()
        QuillenLocalContribution[]
    end
    original_input = verification isa QuillenLocalFactorSequenceVerification ?
        verification.original_input :
        local_correction
    replay_metadata = verification isa QuillenLocalFactorSequenceVerification ?
        verification.replay_metadata :
        (;)
    return QuillenLocalFactorSequenceCertificate(
        original_input,
        R,
        n,
        selected_variable,
        normalized_factors,
        collect(raw_denominators),
        product_denominator,
        local_product,
        local_correction,
        normalized_local_contributions,
        collect(normalized_global_elementary_factors),
        nothing,
        nothing,
        (;),
        replay_metadata,
        verification,
    )
end

function quillen_local_factor_sequence_certificate(
    original_input,
    selected_variable;
    factors,
    local_correction = nothing,
    patched_substitution_witness = nothing,
    chain_witness = nothing,
    witness_metadata = nothing,
    local_evidence = nothing,
    provenance = (;),
    ring = nothing,
    size = nothing,
)
    R, n = _quillen_local_input_ring_size(original_input; ring, size)
    selected = _require_substitution_generator(R, selected_variable)
    raw_factors = collect(factors)
    isempty(raw_factors) &&
        throw(ArgumentError("Quillen local factor sequence certificate requires at least one factor"))
    normalized_factors = [
        _quillen_local_sequence_factor(raw_factor, R, n, index)
        for (index, raw_factor) in enumerate(raw_factors)
    ]
    normalized_local_contributions = _quillen_local_sequence_contributions(normalized_factors, R, n)
    normalized_global_elementary_factors = _quillen_factors(R, n, normalized_local_contributions)
    local_product = _quillen_product(R, n, normalized_global_elementary_factors)
    recorded_local_correction = if local_correction !== nothing
        _quillen_local_require_factor_matrix(local_correction, R, n, "local correction")
    elseif local_evidence !== nothing && hasproperty(local_evidence, :expected_product)
        _quillen_local_require_factor_matrix(local_evidence.expected_product, R, n, "local correction")
    else
        local_product
    end
    stored_witness_metadata = witness_metadata === nothing ?
        (; provenance = provenance, local_evidence = local_evidence) :
        witness_metadata
    raw_denominators = [factor.denominator for factor in normalized_factors]
    product_denominator = prod(raw_denominators; init = one(R))
    placeholder = QuillenLocalFactorSequenceVerification(
        _quillen_local_sequence_original_input(original_input, R, n),
        selected,
        length(normalized_factors),
        raw_denominators,
        product_denominator,
        normalized_local_contributions,
        normalized_global_elementary_factors,
        local_product,
        recorded_local_correction,
        _quillen_denominator_data(normalized_local_contributions),
        _quillen_local_sequence_factor_provenance(normalized_factors),
        false,
        false,
        false,
        false,
        false,
        false,
        (; present = false, ok = false),
        false,
        (;),
        false,
        false,
    )
    provisional = QuillenLocalFactorSequenceCertificate(
        _quillen_local_sequence_original_input(original_input, R, n),
        R,
        n,
        selected,
        normalized_factors,
        raw_denominators,
        product_denominator,
        local_product,
        recorded_local_correction,
        normalized_local_contributions,
        normalized_global_elementary_factors,
        patched_substitution_witness,
        chain_witness,
        stored_witness_metadata,
        nothing,
        placeholder,
    )
    replay_metadata = _quillen_local_sequence_replay_metadata(
        provisional,
        _quillen_denominator_data(normalized_local_contributions),
    )
    provisional = QuillenLocalFactorSequenceCertificate(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.factors,
        provisional.raw_denominators,
        provisional.product_denominator,
        provisional.local_product,
        provisional.local_correction,
        provisional.normalized_local_contributions,
        provisional.normalized_global_elementary_factors,
        provisional.patched_substitution_witness,
        provisional.chain_witness,
        provisional.witness_metadata,
        replay_metadata,
        placeholder,
    )
    verification = replay_quillen_local_factor_sequence(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen local factor sequence certificate data does not replay"))
    return QuillenLocalFactorSequenceCertificate(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.factors,
        provisional.raw_denominators,
        provisional.product_denominator,
        provisional.local_product,
        provisional.local_correction,
        provisional.normalized_local_contributions,
        provisional.normalized_global_elementary_factors,
        provisional.patched_substitution_witness,
        provisional.chain_witness,
        provisional.witness_metadata,
        provisional.replay_metadata,
        verification,
    )
end

function quillen_local_factor_sequence_certificate(
    certificate::QuillenLocalRealizationCertificate;
    factor_provenance = (;
        source = :quillen_local_realization_certificate,
        local_witness_metadata = certificate.witness_metadata,
    ),
    metadata = certificate.witness_metadata,
    chain_witness = nothing,
    local_evidence = nothing,
    provenance = (; source = :length_one_local_factor_sequence),
)
    verify_quillen_local_certificate(certificate) ||
        throw(ArgumentError("Quillen local realization certificate does not replay"))
    factor = QuillenLocalElementaryFactor(
        certificate.correction.row,
        certificate.correction.col,
        certificate.correction.entry,
        certificate.denominator,
        certificate.coverage_multiplier,
        factor_provenance,
        certificate.local_certificate,
        metadata,
    )
    return quillen_local_factor_sequence_certificate(
        certificate.original_input,
        certificate.selected_variable;
        factors = [factor],
        local_correction = certificate.local_correction,
        patched_substitution_witness = certificate.patched_substitution_witness,
        chain_witness = chain_witness,
        witness_metadata = certificate.witness_metadata,
        local_evidence = local_evidence,
        provenance = provenance,
        ring = certificate.ring,
        size = certificate.size,
    )
end

function verify_quillen_local_factor_sequence_certificate(certificate)::Bool
    try
        replay = replay_quillen_local_factor_sequence(certificate)
        return replay.overall_ok &&
               _same_quillen_local_factor_sequence_verification(
                   certificate.verification,
                   replay,
               )
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _murthy_quillen_local_replay(certificate::SL3LocalRealizationCertificate)
    if certificate.branch == :murthy_q0_unit &&
            hasproperty(certificate.witness, :reduction) &&
            certificate.witness.reduction isa SL3LocalMurthyQUnitLocalReduction
        return certificate.witness.reduction.local_factor_replay
    elseif certificate.branch == :murthy_q0_nonunit_bezout_resultant &&
            hasproperty(certificate.witness, :reduction) &&
            hasproperty(certificate.witness.reduction, :local_factor_replay) &&
            certificate.witness.reduction.local_factor_replay !== nothing
        return certificate.witness.reduction.local_factor_replay
    elseif all(factor -> factor isa SL3LocalElementaryFactor, certificate.factors)
        return sl3_local_elementary_factor_replay(
            certificate.target,
            SL3LocalElementaryFactor[certificate.factors...],
            certificate.selected_variable,
        )
    else
        records = map(collect(certificate.factors)) do factor
            factor == identity_matrix(base_ring(factor), nrows(factor)) ?
                sl3_local_elementary_factor(
                    1,
                    2,
                    zero(base_ring(factor)),
                    one(base_ring(factor)),
                    certificate.selected_variable,
                ) :
                only(sl3_local_denominator_one_records_from_matrices([factor], certificate.selected_variable))
        end
        return sl3_local_elementary_factor_replay(
            certificate.target,
            records,
            certificate.selected_variable,
        )
    end
end

function _murthy_quillen_local_record_certificate(record::SL3LocalElementaryFactor)
    return LocalCertificate([record.row, record.col], [record.denominator, record.denominator])
end

function _murthy_quillen_local_sequence_factor(
        record::SL3LocalElementaryFactor,
        index::Int,
        witness_metadata,
)
    return QuillenLocalElementaryFactor(
        record.row,
        record.col,
        record.numerator,
        record.denominator,
        one(record.R),
        _murthy_quillen_local_record_certificate(record),
        (;
            source = :murthy_local_sl3,
            factor_index = index,
            murthy_denominator = record.denominator,
            murthy_local_unit_witness = record.local_unit_witness,
        ),
        (;
            source = :murthy_quillen_local_adapter,
            witness_metadata,
        ),
    )
end

function _murthy_quillen_local_replay_metadata(certificate, replay, mode, witness_metadata)
    return (;
        source = :murthy_quillen_local_adapter,
        murthy_branch = certificate.branch,
        replay_mode = replay.mode,
        adapter_mode = mode,
        selected_variable = replay.selected_variable,
        factor_count = length(replay.factors),
        denominator_product = replay.denominator_product,
        cleared_product = replay.cleared_product,
        ordinary_materialized = replay.materialized_factors !== nothing,
        witness_metadata,
    )
end

function _murthy_quillen_local_factor_sequence(
        original_input,
        selected_variable,
        replay::SL3LocalElementaryFactorReplay,
        witness_metadata,
        provenance,
)
    replay.mode == :ordinary ||
        throw(ArgumentError("Murthy local replay is not materializable over the ordinary base ring"))
    factors = [
        _murthy_quillen_local_sequence_factor(record, index, witness_metadata)
        for (index, record) in enumerate(replay.factors)
    ]
    return quillen_local_factor_sequence_certificate(
        original_input,
        selected_variable;
        factors,
        local_correction = replay.target,
        witness_metadata,
        local_evidence = (;
            source = :murthy_local_sl3,
            replay_mode = replay.mode,
            denominator_product = replay.denominator_product,
            expected_product = replay.target,
        ),
        provenance,
    )
end

function _murthy_quillen_local_single_realization(
        original_input,
        selected_variable,
        replay::SL3LocalElementaryFactorReplay,
        witness_metadata,
)
    replay.mode == :ordinary || return nothing
    record = if length(replay.factors) == 1
        only(replay.factors)
    else
        nonzero_records = filter(record -> record.numerator != zero(record.R), replay.factors)
        length(nonzero_records) == 1 || return nothing
        only(nonzero_records)
    end
    return quillen_local_realization_certificate(
        original_input,
        selected_variable;
        local_certificate = _murthy_quillen_local_record_certificate(record),
        denominator = record.denominator,
        coverage_multiplier = one(record.R),
        correction = QuillenElementaryCorrection(record.row, record.col, record.numerator),
        factors = replay.materialized_factors,
        local_correction = replay.target,
        witness_metadata,
    )
end

function _murthy_quillen_local_adapter(
        certificate::SL3LocalRealizationCertificate,
        original_input,
        selected_variable;
        witness_metadata = (;),
        provenance = (; source = :murthy_quillen_local_adapter),
)
    verify_sl3_local_realization(certificate) ||
        throw(ArgumentError("Murthy local certificate does not replay"))
    R, n = _quillen_local_input_ring_size(original_input)
    n == 3 || throw(ArgumentError("Murthy local Quillen adapter requires a 3x3 input"))
    _same_base_ring(R, base_ring(certificate.target)) ||
        throw(ArgumentError("Murthy local Quillen adapter ring mismatch"))
    original_input == certificate.target ||
        throw(ArgumentError("Murthy local Quillen adapter requires the original input to match the Murthy target"))
    selected = _require_substitution_generator(R, selected_variable)
    selected == certificate.selected_variable ||
        throw(ArgumentError("Murthy local Quillen adapter selected variable mismatch"))

    replay = _murthy_quillen_local_replay(certificate)
    verify_sl3_local_elementary_factor_replay(replay) ||
        throw(ArgumentError("Murthy local factor replay does not verify"))
    replay.target == certificate.target ||
        throw(ArgumentError("Murthy local replay target mismatch"))
    replay.selected_variable == selected ||
        throw(ArgumentError("Murthy local replay selected variable mismatch"))
    (
        replay.factors == certificate.factors ||
        replay.materialized_factors == certificate.factors
    ) ||
        throw(ArgumentError("Murthy local replay factor mismatch"))

    mode = replay.mode == :ordinary ?
        :ordinary_quillen_factor_sequence :
        :localized_replay_handoff
    quillen_sequence = mode == :ordinary_quillen_factor_sequence ?
        _murthy_quillen_local_factor_sequence(
            original_input,
            selected,
            replay,
            witness_metadata,
            provenance,
        ) :
        nothing
    quillen_local = mode == :ordinary_quillen_factor_sequence ?
        _murthy_quillen_local_single_realization(
            original_input,
            selected,
            replay,
            witness_metadata,
        ) :
        nothing
    local_product = quillen_sequence === nothing ? nothing : quillen_sequence.local_product
    local_correction = quillen_sequence === nothing ? replay.target : quillen_sequence.local_correction
    replay_metadata = _murthy_quillen_local_replay_metadata(
        certificate,
        replay,
        mode,
        witness_metadata,
    )
    provisional = MurthyQuillenLocalAdapter(
        original_input,
        R,
        n,
        selected,
        certificate,
        replay,
        mode,
        replay.materialized_factors,
        local_product,
        local_correction,
        quillen_sequence,
        quillen_local,
        witness_metadata,
        replay_metadata,
        nothing,
    )
    verification = _murthy_quillen_local_adapter_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Murthy local Quillen adapter data does not replay"))
    return MurthyQuillenLocalAdapter(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.murthy_certificate,
        provisional.local_factor_replay,
        provisional.mode,
        provisional.materialized_factors,
        provisional.local_product,
        provisional.local_correction,
        provisional.quillen_factor_sequence,
        provisional.quillen_local_certificate,
        provisional.witness_metadata,
        provisional.replay_metadata,
        verification,
    )
end

function _murthy_quillen_local_adapter_summary(adapter::MurthyQuillenLocalAdapter)
    certificate_ok = verify_sl3_local_realization(adapter.murthy_certificate)
    replay_ok = verify_sl3_local_elementary_factor_replay(adapter.local_factor_replay)
    input_ok = adapter.original_input == adapter.murthy_certificate.target
    selected_variable_ok = adapter.selected_variable == adapter.murthy_certificate.selected_variable
    replay_alignment_ok =
        adapter.local_factor_replay.target == adapter.murthy_certificate.target &&
        adapter.local_factor_replay.selected_variable == adapter.selected_variable &&
        (
            adapter.local_factor_replay.factors == adapter.murthy_certificate.factors ||
            adapter.local_factor_replay.materialized_factors == adapter.murthy_certificate.factors
        )
    expected_mode = adapter.local_factor_replay.mode == :ordinary ?
        :ordinary_quillen_factor_sequence :
        :localized_replay_handoff
    mode_ok = adapter.mode == expected_mode
    materialized_ok = adapter.materialized_factors == adapter.local_factor_replay.materialized_factors
    sequence_ok =
        adapter.mode == :ordinary_quillen_factor_sequence ?
        adapter.quillen_factor_sequence isa QuillenLocalFactorSequenceCertificate &&
            verify_quillen_local_factor_sequence_certificate(adapter.quillen_factor_sequence) &&
            adapter.quillen_factor_sequence.original_input == adapter.original_input &&
            adapter.quillen_factor_sequence.selected_variable == adapter.selected_variable &&
            adapter.quillen_factor_sequence.local_product == adapter.local_factor_replay.target &&
            adapter.quillen_factor_sequence.local_correction == adapter.local_factor_replay.target :
        adapter.quillen_factor_sequence === nothing
    expected_local_certificate =
        adapter.mode == :ordinary_quillen_factor_sequence &&
                adapter.quillen_local_certificate !== nothing ?
            _murthy_quillen_local_single_realization(
                adapter.original_input,
                adapter.selected_variable,
                adapter.local_factor_replay,
                adapter.witness_metadata,
            ) :
            nothing
    local_certificate_ok =
        adapter.mode == :localized_replay_handoff ?
        adapter.quillen_local_certificate === nothing :
        (
            adapter.quillen_local_certificate === nothing ||
            (
                adapter.quillen_local_certificate isa QuillenLocalRealizationCertificate &&
                verify_quillen_local_certificate(adapter.quillen_local_certificate) &&
                expected_local_certificate !== nothing &&
                _same_quillen_local_certificate_data(
                    adapter.quillen_local_certificate,
                    expected_local_certificate,
                )
            )
        )
    product_ok =
        adapter.mode == :ordinary_quillen_factor_sequence ?
        adapter.local_product == adapter.local_factor_replay.target :
        adapter.local_product === nothing
    correction_ok = adapter.local_correction == adapter.local_factor_replay.target
    expected_metadata = _murthy_quillen_local_replay_metadata(
        adapter.murthy_certificate,
        adapter.local_factor_replay,
        adapter.mode,
        adapter.witness_metadata,
    )
    replay_metadata_ok = adapter.replay_metadata == expected_metadata
    overall_ok =
        certificate_ok &&
        replay_ok &&
        input_ok &&
        selected_variable_ok &&
        replay_alignment_ok &&
        mode_ok &&
        materialized_ok &&
        sequence_ok &&
        local_certificate_ok &&
        product_ok &&
        correction_ok &&
        replay_metadata_ok
    return (;
        certificate_ok,
        replay_ok,
        input_ok,
        selected_variable_ok,
        replay_alignment_ok,
        mode_ok,
        materialized_ok,
        sequence_ok,
        local_certificate_ok,
        product_ok,
        correction_ok,
        replay_metadata_ok,
        overall_ok,
    )
end

function _verify_murthy_quillen_local_adapter(adapter)::Bool
    try
        replay = _murthy_quillen_local_adapter_summary(adapter)
        return replay.overall_ok && adapter.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _murthy_quillen_local_factor_sequence_certificate(
        adapter::MurthyQuillenLocalAdapter,
)
    _verify_murthy_quillen_local_adapter(adapter) ||
        throw(ArgumentError("Murthy local Quillen adapter does not replay"))
    adapter.quillen_factor_sequence !== nothing ||
        throw(ArgumentError("Murthy local adapter contains localized denominator-cleared replay; #183 must define the localized Quillen local certificate shape before ordinary factor conversion"))
    return adapter.quillen_factor_sequence
end

function _murthy_quillen_local_realization_certificate(
        adapter::MurthyQuillenLocalAdapter,
)
    _verify_murthy_quillen_local_adapter(adapter) ||
        throw(ArgumentError("Murthy local Quillen adapter does not replay"))
    adapter.quillen_local_certificate !== nothing ||
        throw(ArgumentError("Murthy local adapter does not contain a length-one ordinary Quillen local realization certificate"))
    return adapter.quillen_local_certificate
end

function replay_quillen_denominator_cover_candidate(
    candidate::QuillenDenominatorCoverCandidate,
)
    local_certificates = candidate.local_certificates
    isempty(local_certificates) &&
        throw(ArgumentError("Quillen denominator-cover candidate requires at least one local certificate"))
    R = _require_supported_quillen_ring(candidate.ring)
    n = Int(candidate.size)
    n >= 2 || throw(ArgumentError("candidate size must be at least 2"))
    selected_variable = _require_substitution_generator(R, candidate.selected_variable)
    original_input = _quillen_local_sequence_original_input(candidate.original_input, R, n)

    local_certificates_ok = all(verify_quillen_local_factor_sequence_certificate, local_certificates)
    same_ring_ok = all(cert -> cert.ring == R, local_certificates)
    same_size_ok = all(cert -> cert.size == n, local_certificates)
    same_selected_variable_ok = same_ring_ok && all(
        cert -> cert.selected_variable == selected_variable,
        local_certificates,
    )
    same_original_input_ok = same_ring_ok && same_size_ok && all(
        cert -> cert.original_input == original_input,
        local_certificates,
    )

    local_supports = [
        _quillen_local_denominator_support(certificate, index)
        for (index, certificate) in enumerate(local_certificates)
    ]
    replayed_supports_ok = all(support -> support.replay_ok, local_supports)
    raw_denominators = [support.support_denominator for support in local_supports]
    local_supports_ok =
        replayed_supports_ok &&
        _same_quillen_local_denominator_supports(candidate.local_supports, local_supports)
    raw_denominators_ok = candidate.raw_denominators == raw_denominators
    replay_metadata = _quillen_denominator_cover_candidate_replay_metadata(
        local_certificates,
        local_supports,
        raw_denominators,
    )
    replay_metadata_ok = candidate.replay_metadata == replay_metadata
    overall_ok =
        local_certificates_ok &&
        same_original_input_ok &&
        same_ring_ok &&
        same_size_ok &&
        same_selected_variable_ok &&
        local_supports_ok &&
        raw_denominators_ok &&
        replay_metadata_ok

    return QuillenDenominatorCoverCandidateVerification(
        length(local_certificates),
        raw_denominators,
        local_certificates_ok,
        same_original_input_ok,
        same_ring_ok,
        same_size_ok,
        same_selected_variable_ok,
        local_supports,
        local_supports_ok,
        raw_denominators_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function verify_quillen_denominator_cover_candidate(candidate)::Bool
    try
        replay = replay_quillen_denominator_cover_candidate(candidate)
        return replay.overall_ok &&
               _same_quillen_denominator_cover_candidate_verification(
                   candidate.verification,
                   replay,
               )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && return false
        return false
    end
end

function _quillen_denominator_cover_solver_verification(
    source_candidate,
    R,
    raw_denominators,
    exponent::Int,
    coverage_multipliers,
    cover_certificate,
)
    raw_denominator_snapshot = collect(raw_denominators)
    multiplier_snapshot = collect(coverage_multipliers)
    raw_denominator_count = length(raw_denominator_snapshot)
    multiplier_count = length(multiplier_snapshot)
    exact_ring_ok = Oscar.is_exact_type(typeof(zero(coefficient_ring(R))))
    exponent_ok = exponent >= 1
    raw_parent_ring_ok =
        all(denominator -> parent(denominator) == R, raw_denominator_snapshot)
    multiplier_parent_ring_ok =
        all(multiplier -> parent(multiplier) == R, multiplier_snapshot)
    parent_ring_ok =
        raw_denominator_count == multiplier_count &&
        raw_parent_ring_ok &&
        multiplier_parent_ring_ok
    powered_denominators = raw_parent_ring_ok && exponent_ok ?
        [denominator^exponent for denominator in raw_denominator_snapshot] :
        Any[]
    coverage_terms = parent_ring_ok && exponent_ok ?
        [
            multiplier_snapshot[idx] * powered_denominators[idx]
            for idx in eachindex(powered_denominators)
        ] :
        Any[]
    coverage_sum = parent_ring_ok && exponent_ok ?
        sum(coverage_terms; init = zero(R)) :
        zero(R)
    coverage_ok = parent_ring_ok && exact_ring_ok && exponent_ok && coverage_sum == one(R)
    source_candidate_ok =
        source_candidate === nothing ||
        (
            source_candidate isa QuillenDenominatorCoverCandidate &&
            verify_quillen_denominator_cover_candidate(source_candidate) &&
            source_candidate.ring == R &&
            source_candidate.raw_denominators == raw_denominator_snapshot
        )
    cover_certificate_ok =
        cover_certificate isa QuillenDenominatorCoverCertificate &&
        verify_quillen_denominator_cover(cover_certificate)
    cover_certificate_matches =
        cover_certificate_ok &&
        cover_certificate.ring == R &&
        cover_certificate.denominators == powered_denominators &&
        cover_certificate.coverage_multipliers == multiplier_snapshot &&
        cover_certificate.coverage_sum == coverage_sum
    overall_ok = coverage_ok && source_candidate_ok && cover_certificate_matches
    return QuillenDenominatorCoverSolverVerification(
        raw_denominator_count,
        multiplier_count,
        raw_denominator_snapshot,
        exponent,
        powered_denominators,
        multiplier_snapshot,
        parent_ring_ok,
        exact_ring_ok,
        exponent_ok,
        source_candidate_ok,
        coverage_terms,
        coverage_sum,
        coverage_ok,
        cover_certificate_ok,
        cover_certificate_matches,
        overall_ok,
    )
end

function _same_quillen_denominator_cover_solver_verification(
    left::QuillenDenominatorCoverSolverVerification,
    right::QuillenDenominatorCoverSolverVerification,
)::Bool
    return left.raw_denominator_count == right.raw_denominator_count &&
           left.multiplier_count == right.multiplier_count &&
           left.raw_denominators == right.raw_denominators &&
           left.exponent == right.exponent &&
           left.powered_denominators == right.powered_denominators &&
           left.coverage_multipliers == right.coverage_multipliers &&
           left.parent_ring_ok == right.parent_ring_ok &&
           left.exact_ring_ok == right.exact_ring_ok &&
           left.exponent_ok == right.exponent_ok &&
           left.source_candidate_ok == right.source_candidate_ok &&
           left.coverage_terms == right.coverage_terms &&
           left.coverage_sum == right.coverage_sum &&
           left.coverage_ok == right.coverage_ok &&
           left.cover_certificate_ok == right.cover_certificate_ok &&
           left.cover_certificate_matches == right.cover_certificate_matches &&
           left.overall_ok == right.overall_ok
end

function replay_quillen_denominator_cover_solver_result(
    result::QuillenDenominatorCoverSolverResult,
)
    _require_quillen_denominator_cover_ring(result.ring)
    return _quillen_denominator_cover_solver_verification(
        result.source_candidate,
        result.ring,
        result.raw_denominators,
        result.exponent,
        result.coverage_multipliers,
        result.cover_certificate,
    )
end

function verify_quillen_denominator_cover_solver_result(result)::Bool
    try
        replay = replay_quillen_denominator_cover_solver_result(result)
        return replay.overall_ok &&
               result.raw_denominators == replay.raw_denominators &&
               result.exponent == replay.exponent &&
               result.powered_denominators == replay.powered_denominators &&
               result.coverage_multipliers == replay.coverage_multipliers &&
               result.coverage_terms == replay.coverage_terms &&
               result.coverage_sum == replay.coverage_sum &&
               _same_quillen_denominator_cover_solver_verification(
                   result.verification,
                   replay,
               )
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _quillen_cover_exponent_range(max_exponent::Integer, exponent)
    bound = Int(max_exponent)
    bound >= 1 || throw(ArgumentError("coverage not proven: max_exponent must be positive"))
    if exponent === nothing
        return 1:bound
    end
    exponent isa Integer ||
        throw(ArgumentError("coverage not proven: exponent must be an integer"))
    chosen = Int(exponent)
    1 <= chosen <= bound ||
        throw(ArgumentError("coverage not proven: requested exponent is outside the configured bound"))
    return chosen:chosen
end

function _quillen_supplied_cover_multipliers(coverage_multipliers, supplied_multipliers)
    coverage_multipliers !== nothing && supplied_multipliers !== nothing &&
        throw(ArgumentError("coverage not proven: provide only one supplied multiplier collection"))
    return coverage_multipliers === nothing ? supplied_multipliers : coverage_multipliers
end

function _quillen_cover_coordinate_multipliers(R, powered_denominators)
    cover_ideal = ideal(R, powered_denominators)
    one(R) in cover_ideal || return nothing
    coordinates = Oscar.coordinates(one(R), cover_ideal)
    return [R(coordinates[1, idx]) for idx in eachindex(powered_denominators)]
end

function _quillen_denominator_cover_solver_result(
    source_candidate,
    R,
    raw_denominators,
    exponent::Int,
    coverage_multipliers,
)
    cover_certificate = quillen_denominator_cover_certificate(
        R,
        [denominator^exponent for denominator in raw_denominators],
        coverage_multipliers,
    )
    verification = _quillen_denominator_cover_solver_verification(
        source_candidate,
        R,
        raw_denominators,
        exponent,
        coverage_multipliers,
        cover_certificate,
    )
    verification.overall_ok ||
        throw(ArgumentError("coverage not proven: denominator-cover solver result does not replay"))
    return QuillenDenominatorCoverSolverResult(
        source_candidate,
        R,
        verification.raw_denominators,
        verification.exponent,
        verification.powered_denominators,
        verification.coverage_multipliers,
        verification.coverage_terms,
        verification.coverage_sum,
        cover_certificate,
        verification,
    )
end

function solve_quillen_denominator_cover(
    R,
    raw_denominators;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
    source_candidate = nothing,
)
    _require_quillen_denominator_cover_ring(R)
    normalized_denominators =
        _normalize_quillen_cover_elements(R, raw_denominators, "cover denominator")
    isempty(normalized_denominators) &&
        throw(ArgumentError("coverage not proven: denominator list must be nonempty"))
    exponent_range = _quillen_cover_exponent_range(max_exponent, exponent)
    supplied = _quillen_supplied_cover_multipliers(
        coverage_multipliers,
        supplied_multipliers,
    )
    normalized_supplied = supplied === nothing ?
        nothing :
        _normalize_quillen_cover_elements(R, supplied, "cover coverage multiplier")

    for chosen_exponent in exponent_range
        powered_denominators = [
            denominator^chosen_exponent for denominator in normalized_denominators
        ]
        multipliers = if normalized_supplied === nothing
            _quillen_cover_coordinate_multipliers(R, powered_denominators)
        else
            normalized_supplied
        end
        multipliers === nothing && continue
        try
            result = _quillen_denominator_cover_solver_result(
                source_candidate,
                R,
                normalized_denominators,
                Int(chosen_exponent),
                multipliers,
            )
            verify_quillen_denominator_cover_solver_result(result) && return result
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
        end
    end

    throw(ArgumentError("coverage not proven: bounded denominator-cover search failed"))
end

function solve_quillen_denominator_cover(
    candidate::QuillenDenominatorCoverCandidate;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
)
    verify_quillen_denominator_cover_candidate(candidate) ||
        throw(ArgumentError("coverage not proven: denominator-cover candidate does not replay"))
    return solve_quillen_denominator_cover(
        candidate.ring,
        candidate.raw_denominators;
        max_exponent,
        exponent,
        coverage_multipliers,
        supplied_multipliers,
        source_candidate = candidate,
    )
end

function _quillen_patch_substitution_step_metadata(
    solver_result::QuillenDenominatorCoverSolverResult,
    selected_variable,
    sign_convention::Symbol,
    step_index::Int,
)
    return (;
        source = :park_woodburn_substitution_chain,
        step_index = step_index,
        selected_variable = selected_variable,
        raw_denominator = solver_result.raw_denominators[step_index],
        exponent = solver_result.exponent,
        powered_denominator = solver_result.powered_denominators[step_index],
        coverage_multiplier = solver_result.coverage_multipliers[step_index],
        coverage_term = solver_result.coverage_terms[step_index],
        sign_convention = sign_convention,
    )
end

function _quillen_patch_substitution_chain_metadata(
    solver_result::QuillenDenominatorCoverSolverResult,
    selected_variable,
    sign_convention::Symbol,
    metadata,
)
    return (;
        source = :park_woodburn_substitution_chain,
        selected_variable = selected_variable,
        exponent = solver_result.exponent,
        denominator_count = length(solver_result.raw_denominators),
        coverage_sum = solver_result.coverage_sum,
        sign_convention = sign_convention,
        metadata = metadata,
    )
end

function _quillen_patch_next_coefficient(
    previous,
    powered_denominator,
    coverage_multiplier,
    sign_convention::Symbol,
)
    sign_convention == :park_woodburn_minus ||
        throw(ArgumentError("unsupported Park-Woodburn substitution sign convention"))
    return previous - coverage_multiplier * powered_denominator
end

function _same_quillen_patch_substitution_step(
    left::QuillenPatchSubstitutionStep,
    right::QuillenPatchSubstitutionStep,
)::Bool
    return left.step_index == right.step_index &&
           left.selected_variable == right.selected_variable &&
           left.raw_denominator == right.raw_denominator &&
           left.exponent == right.exponent &&
           left.powered_denominator == right.powered_denominator &&
           left.coverage_multiplier == right.coverage_multiplier &&
           left.sign_convention == right.sign_convention &&
           left.previous_coefficient == right.previous_coefficient &&
           left.next_coefficient == right.next_coefficient &&
           left.previous_matrix == right.previous_matrix &&
           left.next_matrix == right.next_matrix &&
           left.bracket_target == right.bracket_target &&
           left.replay_metadata == right.replay_metadata
end

function _same_quillen_patch_substitution_steps(left, right)::Bool
    length(left) == length(right) || return false
    return all(
        _same_quillen_patch_substitution_step(left[idx], right[idx])
        for idx in eachindex(left)
    )
end

function _quillen_patch_expected_substitution_chain(
    A,
    selected_variable,
    solver_result::QuillenDenominatorCoverSolverResult,
    sign_convention::Symbol,
    metadata,
)
    R = base_ring(A)
    cumulative_coefficients = Any[one(R)]
    intermediate_matrices = Any[
        _quillen_substitute_matrix_scaled_variable(A, selected_variable, one(R)),
    ]
    steps = QuillenPatchSubstitutionStep[]
    bracket_matrices = Any[]

    for idx in eachindex(solver_result.raw_denominators)
        previous_coefficient = last(cumulative_coefficients)
        next_coefficient = _quillen_patch_next_coefficient(
            previous_coefficient,
            solver_result.powered_denominators[idx],
            solver_result.coverage_multipliers[idx],
            sign_convention,
        )
        previous_matrix = last(intermediate_matrices)
        next_matrix = _quillen_substitute_matrix_scaled_variable(
            A,
            selected_variable,
            next_coefficient,
        )
        bracket_target = inv(previous_matrix) * next_matrix
        push!(cumulative_coefficients, next_coefficient)
        push!(intermediate_matrices, next_matrix)
        push!(bracket_matrices, bracket_target)
        push!(
            steps,
            QuillenPatchSubstitutionStep(
                idx,
                selected_variable,
                solver_result.raw_denominators[idx],
                solver_result.exponent,
                solver_result.powered_denominators[idx],
                solver_result.coverage_multipliers[idx],
                sign_convention,
                previous_coefficient,
                next_coefficient,
                previous_matrix,
                next_matrix,
                bracket_target,
                _quillen_patch_substitution_step_metadata(
                    solver_result,
                    selected_variable,
                    sign_convention,
                    idx,
                ),
            ),
        )
    end

    base_term = _quillen_substitute_matrix_scaled_variable(A, selected_variable, zero(R))
    telescope_product = A
    for bracket in bracket_matrices
        telescope_product *= bracket
    end

    return (;
        cumulative_coefficients = cumulative_coefficients,
        intermediate_matrices = intermediate_matrices,
        steps = steps,
        bracket_matrices = bracket_matrices,
        final_coefficient = last(cumulative_coefficients),
        base_term = base_term,
        telescope_product = telescope_product,
        replay_metadata = _quillen_patch_substitution_chain_metadata(
            solver_result,
            selected_variable,
            sign_convention,
            metadata,
        ),
    )
end

function _same_quillen_patch_substitution_chain_verification(
    left::QuillenPatchSubstitutionChainVerification,
    right::QuillenPatchSubstitutionChainVerification,
)::Bool
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
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _quillen_patch_substitution_chain_verification(
    A,
    R,
    size::Int,
    selected_variable,
    sign_convention::Symbol,
    solver_result::QuillenDenominatorCoverSolverResult,
    cumulative_coefficients,
    intermediate_matrices,
    steps,
    bracket_matrices,
    base_term,
    metadata,
    replay_metadata,
)
    solver_result_ok = verify_quillen_denominator_cover_solver_result(solver_result)
    ring_ok = solver_result_ok && solver_result.ring == R
    matrix_ok = nrows(A) == size && ncols(A) == size && base_ring(A) == R
    selected = _require_substitution_generator(R, selected_variable)
    selected_variable_ok = selected == selected_variable
    sign_convention_ok = sign_convention == :park_woodburn_minus
    expected = _quillen_patch_expected_substitution_chain(
        A,
        selected,
        solver_result,
        sign_convention,
        metadata,
    )

    coefficient_count = length(cumulative_coefficients)
    step_count = length(steps)
    bracket_count = length(bracket_matrices)
    cumulative_coefficients_ok =
        collect(cumulative_coefficients) == expected.cumulative_coefficients
    intermediate_matrices_ok =
        collect(intermediate_matrices) == expected.intermediate_matrices
    steps_ok = _same_quillen_patch_substitution_steps(steps, expected.steps)
    bracket_matrices_ok = collect(bracket_matrices) == expected.bracket_matrices
    final_coefficient_ok = expected.final_coefficient == zero(R)
    base_term_ok = base_term == expected.base_term &&
                   base_term == last(expected.intermediate_matrices)
    telescope_ok = expected.telescope_product == expected.base_term
    replay_metadata_ok = replay_metadata == expected.replay_metadata
    overall_ok =
        solver_result_ok &&
        ring_ok &&
        matrix_ok &&
        selected_variable_ok &&
        sign_convention_ok &&
        coefficient_count == length(expected.cumulative_coefficients) &&
        step_count == length(expected.steps) &&
        bracket_count == length(expected.bracket_matrices) &&
        cumulative_coefficients_ok &&
        intermediate_matrices_ok &&
        steps_ok &&
        bracket_matrices_ok &&
        final_coefficient_ok &&
        base_term_ok &&
        telescope_ok &&
        replay_metadata_ok

    return QuillenPatchSubstitutionChainVerification(
        solver_result_ok,
        ring_ok,
        matrix_ok,
        selected_variable_ok,
        sign_convention_ok,
        coefficient_count,
        step_count,
        bracket_count,
        expected.cumulative_coefficients,
        cumulative_coefficients_ok,
        expected.intermediate_matrices,
        intermediate_matrices_ok,
        expected.steps,
        steps_ok,
        expected.bracket_matrices,
        bracket_matrices_ok,
        expected.final_coefficient,
        final_coefficient_ok,
        expected.base_term,
        base_term_ok,
        expected.telescope_product,
        telescope_ok,
        expected.replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function replay_quillen_patch_substitution_chain(
    chain::QuillenPatchSubstitutionChain,
)
    R = _require_quillen_denominator_cover_ring(chain.ring)
    _require_square_matrix(chain.original_matrix, "substitution-chain original matrix") ==
        chain.size || throw(DimensionMismatch("substitution-chain size must match original matrix"))
    return _quillen_patch_substitution_chain_verification(
        chain.original_matrix,
        R,
        chain.size,
        chain.selected_variable,
        chain.sign_convention,
        chain.solver_result,
        chain.cumulative_coefficients,
        chain.intermediate_matrices,
        chain.steps,
        chain.bracket_matrices,
        chain.base_term,
        chain.metadata,
        chain.replay_metadata,
    )
end

function verify_quillen_patch_substitution_chain(chain)::Bool
    try
        replay = replay_quillen_patch_substitution_chain(chain)
        return replay.overall_ok &&
               _same_quillen_patch_substitution_chain_verification(
                   chain.verification,
                   replay,
               )
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function quillen_patch_substitution_chain(
    A,
    selected_variable,
    solver_result::QuillenDenominatorCoverSolverResult;
    sign_convention::Symbol = :park_woodburn_minus,
    metadata = (;),
)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    solver_result.ring == R ||
        throw(ArgumentError("substitution-chain solver result ring must match matrix ring"))
    verify_quillen_denominator_cover_solver_result(solver_result) ||
        throw(ArgumentError("substitution-chain solver result must replay"))
    n = _require_square_matrix(A, "substitution-chain original matrix")
    selected = _require_substitution_generator(R, selected_variable)
    sign_convention == :park_woodburn_minus ||
        throw(ArgumentError("unsupported Park-Woodburn substitution sign convention"))

    expected = _quillen_patch_expected_substitution_chain(
        A,
        selected,
        solver_result,
        sign_convention,
        metadata,
    )
    verification = _quillen_patch_substitution_chain_verification(
        A,
        R,
        n,
        selected,
        sign_convention,
        solver_result,
        expected.cumulative_coefficients,
        expected.intermediate_matrices,
        expected.steps,
        expected.bracket_matrices,
        expected.base_term,
        metadata,
        expected.replay_metadata,
    )
    verification.overall_ok ||
        throw(ArgumentError("Park-Woodburn substitution chain does not replay"))
    return QuillenPatchSubstitutionChain(
        A,
        R,
        n,
        selected,
        sign_convention,
        solver_result,
        expected.cumulative_coefficients,
        expected.intermediate_matrices,
        expected.steps,
        expected.bracket_matrices,
        expected.base_term,
        metadata,
        expected.replay_metadata,
        verification,
    )
end

function extract_quillen_denominator_cover_candidate(certificates)
    local_certificates = collect(certificates)
    isempty(local_certificates) &&
        throw(ArgumentError("Quillen denominator-cover candidate requires at least one local certificate"))
    all(verify_quillen_local_factor_sequence_certificate, local_certificates) ||
        throw(ArgumentError("Quillen denominator-cover candidate requires verified local certificates"))

    first_certificate = first(local_certificates)
    R = _require_supported_quillen_ring(first_certificate.ring)
    n = Int(first_certificate.size)
    n >= 2 || throw(ArgumentError("candidate size must be at least 2"))
    selected_variable = _require_substitution_generator(R, first_certificate.selected_variable)
    original_input = _quillen_local_sequence_original_input(first_certificate.original_input, R, n)

    all(cert -> cert.ring == R, local_certificates) ||
        throw(ArgumentError("Quillen denominator-cover candidate requires matching rings"))
    all(cert -> cert.size == n, local_certificates) ||
        throw(ArgumentError("Quillen denominator-cover candidate requires matching sizes"))
    all(cert -> cert.selected_variable == selected_variable, local_certificates) ||
        throw(ArgumentError("Quillen denominator-cover candidate requires matching selected variables"))
    all(cert -> cert.original_input == original_input, local_certificates) ||
        throw(ArgumentError("Quillen denominator-cover candidate requires matching original inputs"))

    local_supports = [
        _quillen_local_denominator_support(certificate, index)
        for (index, certificate) in enumerate(local_certificates)
    ]
    raw_denominators = [support.support_denominator for support in local_supports]
    replay_metadata = _quillen_denominator_cover_candidate_replay_metadata(
        local_certificates,
        local_supports,
        raw_denominators,
    )
    placeholder = QuillenDenominatorCoverCandidateVerification(
        length(local_certificates),
        raw_denominators,
        false,
        false,
        false,
        false,
        false,
        local_supports,
        false,
        false,
        replay_metadata,
        false,
        false,
    )
    provisional = QuillenDenominatorCoverCandidate(
        original_input,
        R,
        n,
        selected_variable,
        local_certificates,
        raw_denominators,
        local_supports,
        replay_metadata,
        placeholder,
    )
    verification = replay_quillen_denominator_cover_candidate(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen denominator-cover candidate data does not replay"))
    return QuillenDenominatorCoverCandidate(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.local_certificates,
        provisional.raw_denominators,
        provisional.local_supports,
        provisional.replay_metadata,
        verification,
    )
end

function replay_quillen_local_contribution_normalization(
    normalized::QuillenLocalContributionNormalization,
)
    certificate = normalized.local_certificate
    cover = normalized.cover_certificate
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    selected_variable = _require_substitution_generator(R, normalized.selected_variable)
    local_replay = _quillen_local_certificate_replay_summary(certificate)

    local_certificate_ok = verify_quillen_local_certificate(certificate)
    cover_certificate_ok = verify_quillen_denominator_cover(cover)
    original_input_ok =
        normalized.original_input == certificate.original_input &&
        local_replay.original_input == normalized.original_input
    selected_variable_ok =
        selected_variable == certificate.selected_variable &&
        local_replay.selected_variable == selected_variable
    cover_ring_ok = cover.ring == R
    cover_index_ok = 1 <= normalized.cover_index <= length(cover.denominators)
    cover_pair = cover_index_ok ?
        QuillenDenominatorData(
            cover.denominators[normalized.cover_index],
            cover.coverage_multipliers[normalized.cover_index],
        ) :
        nothing
    cover_pair_ok =
        cover_ring_ok &&
        cover_index_ok &&
        cover_pair.denominator == normalized.denominator &&
        cover_pair.coverage_multiplier == normalized.coverage_multiplier &&
        certificate.denominator == normalized.denominator &&
        certificate.coverage_multiplier == normalized.coverage_multiplier

    patched_substitution = _quillen_local_patched_witness_summary(
        normalized.patched_substitution_witness,
        R,
        n,
        selected_variable,
    )
    patched_substitution_ok =
        normalized.patched_substitution_witness == certificate.patched_substitution_witness &&
        patched_substitution.ok &&
        normalized.patched_substitution == patched_substitution

    replayed_contribution = _normalize_quillen_contribution(
        QuillenLocalContribution(
            certificate.local_certificate,
            normalized.denominator,
            normalized.coverage_multiplier,
            certificate.correction,
        ),
        R,
        n,
    )
    local_contribution_ok = _same_quillen_local_contribution(
        normalized.local_contribution,
        replayed_contribution,
    )
    weighted_global_elementary_factors = _quillen_factors(R, n, [replayed_contribution])
    weighted_global_elementary_factors_ok = _same_quillen_factors(
        normalized.weighted_global_elementary_factors,
        weighted_global_elementary_factors,
    )
    local_product_ok = normalized.local_product == local_replay.local_product
    local_correction_ok = normalized.local_correction == local_replay.local_correction
    replay_metadata = _quillen_normalization_replay_metadata(
        certificate,
        cover,
        normalized.cover_index,
    )
    replay_metadata_ok = normalized.replay_metadata == replay_metadata
    overall_ok =
        local_certificate_ok &&
        cover_certificate_ok &&
        original_input_ok &&
        selected_variable_ok &&
        cover_ring_ok &&
        cover_pair_ok &&
        patched_substitution_ok &&
        local_product_ok &&
        local_correction_ok &&
        local_contribution_ok &&
        weighted_global_elementary_factors_ok &&
        replay_metadata_ok

    return QuillenLocalContributionNormalizationVerification(
        local_certificate_ok,
        cover_certificate_ok,
        original_input_ok,
        selected_variable_ok,
        cover_ring_ok,
        cover_index_ok,
        cover_pair,
        cover_pair_ok,
        patched_substitution,
        patched_substitution_ok,
        local_replay.local_product,
        local_product_ok,
        local_replay.local_correction,
        local_correction_ok,
        local_contribution_ok,
        weighted_global_elementary_factors,
        weighted_global_elementary_factors_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function normalize_quillen_local_contribution(
    certificate::QuillenLocalRealizationCertificate,
    cover::QuillenDenominatorCoverCertificate;
    original_input = certificate.original_input,
    selected_variable = certificate.selected_variable,
    patched_substitution_witness = certificate.patched_substitution_witness,
)
    verify_quillen_local_certificate(certificate) ||
        throw(ArgumentError("Quillen local realization certificate does not replay"))
    verify_quillen_denominator_cover(cover) ||
        throw(ArgumentError("Quillen denominator cover certificate does not replay"))
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    cover.ring == R ||
        throw(ArgumentError("Quillen local contribution normalization requires cover and local certificate rings to match"))
    selected = _require_substitution_generator(R, selected_variable)
    selected == certificate.selected_variable ||
        throw(ArgumentError("Quillen local contribution normalization requires the selected variable recorded by the local certificate"))
    original_input == certificate.original_input ||
        throw(ArgumentError("Quillen local contribution normalization requires the original input recorded by the local certificate"))
    patched_substitution_witness == certificate.patched_substitution_witness ||
        throw(ArgumentError("Quillen local contribution normalization requires the patched-substitution witness recorded by the local certificate"))

    cover_index = _quillen_cover_pair_index(
        cover,
        certificate.denominator,
        certificate.coverage_multiplier,
    )
    local_contribution = _normalize_quillen_contribution(
        QuillenLocalContribution(
            certificate.local_certificate,
            certificate.denominator,
            certificate.coverage_multiplier,
            certificate.correction,
        ),
        R,
        n,
    )
    weighted_global_elementary_factors = _quillen_factors(R, n, [local_contribution])
    patched_substitution = _quillen_local_patched_witness_summary(
        patched_substitution_witness,
        R,
        n,
        selected,
    )
    patched_substitution.ok ||
        throw(ArgumentError("Quillen local contribution normalization patched substitution does not replay"))
    replay_metadata = _quillen_normalization_replay_metadata(certificate, cover, cover_index)
    provisional = QuillenLocalContributionNormalization(
        certificate,
        cover,
        original_input,
        selected,
        certificate.denominator,
        certificate.coverage_multiplier,
        cover_index,
        patched_substitution_witness,
        patched_substitution,
        certificate.local_product,
        certificate.local_correction,
        local_contribution,
        weighted_global_elementary_factors,
        replay_metadata,
        QuillenLocalContributionNormalizationVerification(
            false,
            false,
            false,
            false,
            false,
            false,
            nothing,
            false,
            (; present = false, ok = false),
            false,
            certificate.local_product,
            false,
            certificate.local_correction,
            false,
            false,
            Any[],
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_quillen_local_contribution_normalization(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen local contribution normalization data does not replay"))
    return QuillenLocalContributionNormalization(
        provisional.local_certificate,
        provisional.cover_certificate,
        provisional.original_input,
        provisional.selected_variable,
        provisional.denominator,
        provisional.coverage_multiplier,
        provisional.cover_index,
        provisional.patched_substitution_witness,
        provisional.patched_substitution,
        provisional.local_product,
        provisional.local_correction,
        provisional.local_contribution,
        provisional.weighted_global_elementary_factors,
        provisional.replay_metadata,
        verification,
    )
end

function normalize_quillen_local_contributions(
    certificates,
    cover::QuillenDenominatorCoverCertificate;
    original_input = nothing,
    selected_variable = nothing,
)
    return [
        normalize_quillen_local_contribution(
            certificate,
            cover;
            original_input = original_input === nothing ? certificate.original_input : original_input,
            selected_variable = selected_variable === nothing ? certificate.selected_variable : selected_variable,
        )
        for certificate in collect(certificates)
    ]
end

function verify_quillen_local_contribution_normalization(normalized)::Bool
    try
        replay = replay_quillen_local_contribution_normalization(normalized)
        return replay.overall_ok &&
               _same_quillen_normalization_verification(normalized.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

struct QuillenSequenceContributionExpansionVerification
    local_certificate_ok::Bool
    solver_result_ok::Bool
    local_index_ok::Bool
    solver_context_ok::Bool
    powered_denominator
    coverage_multiplier
    cover_term
    cover_term_ok::Bool
    factor_provenance::Vector
    factor_provenance_ok::Bool
    global_elementary_factors::Vector
    global_elementary_factors_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenSequenceContributionExpansion
    local_certificate::QuillenLocalFactorSequenceCertificate
    solver_result::QuillenDenominatorCoverSolverResult
    local_index::Int
    powered_denominator
    coverage_multiplier
    cover_term
    factor_provenance::Vector
    global_elementary_factors::Vector
    replay_metadata
    verification::QuillenSequenceContributionExpansionVerification
end

function _quillen_sequence_expansion_metadata(
    certificate::QuillenLocalFactorSequenceCertificate,
    solver_result::QuillenDenominatorCoverSolverResult,
    local_index::Int,
    factor_provenance,
)
    return (;
        source = :quillen_supplied_local_sequence_expansion,
        local_index = local_index,
        factor_count = length(certificate.factors),
        raw_denominator = solver_result.raw_denominators[local_index],
        powered_denominator = solver_result.powered_denominators[local_index],
        coverage_multiplier = solver_result.coverage_multipliers[local_index],
        cover_term = solver_result.coverage_terms[local_index],
        factor_provenance = factor_provenance,
        local_replay_metadata = certificate.replay_metadata,
    )
end

function _quillen_sequence_expansion_factors(
    certificate::QuillenLocalFactorSequenceCertificate,
    cover_term,
)
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    factor_type = typeof(identity_matrix(R, n))
    factors = factor_type[]
    for factor in certificate.factors
        push!(
            factors,
            elementary_matrix(
                n,
                factor.row,
                factor.col,
                _coerce_into_ring(R, cover_term * factor.numerator, "sequence expansion entry"),
                R,
            ),
        )
    end
    return factors
end

function _quillen_solver_source_candidate_matches_sequence(
    solver_result::QuillenDenominatorCoverSolverResult,
    certificate::QuillenLocalFactorSequenceCertificate,
    local_index::Int,
)::Bool
    source_candidate = solver_result.source_candidate
    source_candidate isa QuillenDenominatorCoverCandidate || return false
    verify_quillen_denominator_cover_candidate(source_candidate) || return false
    1 <= local_index <= length(source_candidate.local_certificates) || return false
    return source_candidate.original_input == certificate.original_input &&
           source_candidate.ring == certificate.ring &&
           source_candidate.size == certificate.size &&
           source_candidate.selected_variable == certificate.selected_variable &&
           source_candidate.local_certificates[local_index] == certificate &&
           source_candidate.raw_denominators == solver_result.raw_denominators
end

function _same_quillen_sequence_expansion_verification(
    left::QuillenSequenceContributionExpansionVerification,
    right::QuillenSequenceContributionExpansionVerification,
)::Bool
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
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function replay_quillen_sequence_contribution_expansion(
    expansion::QuillenSequenceContributionExpansion,
)
    certificate = expansion.local_certificate
    solver_result = expansion.solver_result
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    local_certificate_ok = verify_quillen_local_factor_sequence_certificate(certificate)
    solver_result_ok = verify_quillen_denominator_cover_solver_result(solver_result)
    local_index_ok = 1 <= expansion.local_index <= length(solver_result.raw_denominators)
    solver_context_ok =
        solver_result_ok &&
        solver_result.ring == R &&
        local_index_ok &&
        solver_result.raw_denominators[expansion.local_index] == certificate.product_denominator &&
        _quillen_solver_source_candidate_matches_sequence(
            solver_result,
            certificate,
            expansion.local_index,
        )
    powered_denominator = local_index_ok ? solver_result.powered_denominators[expansion.local_index] : zero(R)
    coverage_multiplier = local_index_ok ? solver_result.coverage_multipliers[expansion.local_index] : zero(R)
    cover_term = local_index_ok ? solver_result.coverage_terms[expansion.local_index] : zero(R)
    cover_term_ok =
        expansion.powered_denominator == powered_denominator &&
        expansion.coverage_multiplier == coverage_multiplier &&
        expansion.cover_term == cover_term &&
        cover_term == coverage_multiplier * powered_denominator
    factor_provenance = _quillen_local_sequence_factor_provenance(certificate.factors)
    factor_provenance_ok = expansion.factor_provenance == factor_provenance
    global_elementary_factors = _quillen_sequence_expansion_factors(certificate, cover_term)
    global_elementary_factors_ok =
        _same_quillen_factors(expansion.global_elementary_factors, global_elementary_factors)
    replay_metadata = _quillen_sequence_expansion_metadata(
        certificate,
        solver_result,
        expansion.local_index,
        factor_provenance,
    )
    replay_metadata_ok = expansion.replay_metadata == replay_metadata
    overall_ok =
        local_certificate_ok &&
        solver_result_ok &&
        local_index_ok &&
        solver_context_ok &&
        cover_term_ok &&
        factor_provenance_ok &&
        global_elementary_factors_ok &&
        replay_metadata_ok
    return QuillenSequenceContributionExpansionVerification(
        local_certificate_ok,
        solver_result_ok,
        local_index_ok,
        solver_context_ok,
        powered_denominator,
        coverage_multiplier,
        cover_term,
        cover_term_ok,
        factor_provenance,
        factor_provenance_ok,
        global_elementary_factors,
        global_elementary_factors_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function quillen_sequence_contribution_expansion(
    certificate::QuillenLocalFactorSequenceCertificate,
    solver_result::QuillenDenominatorCoverSolverResult,
    local_index::Int,
)
    verify_quillen_local_factor_sequence_certificate(certificate) ||
        throw(ArgumentError("Quillen supplied evidence assembly requires verified local sequence certificates"))
    verify_quillen_denominator_cover_solver_result(solver_result) ||
        throw(ArgumentError("Quillen supplied evidence assembly requires a verified denominator-cover solver result"))
    1 <= local_index <= length(solver_result.raw_denominators) ||
        throw(ArgumentError("Quillen supplied evidence assembly local index is outside the solver result"))
    solver_result.ring == certificate.ring &&
        solver_result.raw_denominators[local_index] == certificate.product_denominator ||
        throw(ArgumentError("Quillen supplied evidence assembly solver denominator must match local sequence provenance"))
    factor_provenance = _quillen_local_sequence_factor_provenance(certificate.factors)
    powered_denominator = solver_result.powered_denominators[local_index]
    coverage_multiplier = solver_result.coverage_multipliers[local_index]
    cover_term = solver_result.coverage_terms[local_index]
    global_elementary_factors = _quillen_sequence_expansion_factors(certificate, cover_term)
    replay_metadata = _quillen_sequence_expansion_metadata(
        certificate,
        solver_result,
        local_index,
        factor_provenance,
    )
    provisional = QuillenSequenceContributionExpansion(
        certificate,
        solver_result,
        Int(local_index),
        powered_denominator,
        coverage_multiplier,
        cover_term,
        factor_provenance,
        global_elementary_factors,
        replay_metadata,
        QuillenSequenceContributionExpansionVerification(
            false,
            false,
            false,
            false,
            powered_denominator,
            coverage_multiplier,
            cover_term,
            false,
            factor_provenance,
            false,
            typeof(identity_matrix(certificate.ring, certificate.size))[],
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_quillen_sequence_contribution_expansion(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen supplied evidence sequence expansion does not replay"))
    return QuillenSequenceContributionExpansion(
        provisional.local_certificate,
        provisional.solver_result,
        provisional.local_index,
        provisional.powered_denominator,
        provisional.coverage_multiplier,
        provisional.cover_term,
        provisional.factor_provenance,
        provisional.global_elementary_factors,
        provisional.replay_metadata,
        verification,
    )
end

function verify_quillen_sequence_contribution_expansion(expansion)::Bool
    try
        replay = replay_quillen_sequence_contribution_expansion(expansion)
        return replay.overall_ok &&
               _same_quillen_sequence_expansion_verification(expansion.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

struct QuillenSuppliedEvidencePatchAssemblyVerification
    local_certificates_ok::Bool
    denominator_candidate_ok::Bool
    denominator_candidate_matches::Bool
    solver_result_ok::Bool
    solver_source_candidate_ok::Bool
    cover_certificate_ok::Bool
    substitution_chain_ok::Bool
    substitution_chain_matches::Bool
    base_term_ok::Bool
    sequence_expansions_ok::Bool
    global_elementary_factors::Vector
    global_elementary_factors_ok::Bool
    product
    product_ok::Bool
    target
    target_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenSuppliedEvidencePatchAssembly
    ring
    size::Int
    substitution_variable
    original_input
    local_certificates::Vector{QuillenLocalFactorSequenceCertificate}
    denominator_candidate::QuillenDenominatorCoverCandidate
    solver_result::QuillenDenominatorCoverSolverResult
    cover_certificate::QuillenDenominatorCoverCertificate
    substitution_chain::QuillenPatchSubstitutionChain
    base_term_policy::Symbol
    base_term
    base_term_factors::Vector
    base_term_product
    sequence_expansions::Vector{QuillenSequenceContributionExpansion}
    sequence_elementary_factors::Vector
    global_elementary_factors::Vector
    product
    target
    replay_metadata
    verification::QuillenSuppliedEvidencePatchAssemblyVerification
end

function _quillen_supplied_patch_metadata(
    candidate::QuillenDenominatorCoverCandidate,
    solver_result::QuillenDenominatorCoverSolverResult,
    substitution_chain::QuillenPatchSubstitutionChain,
    base_term_policy::Symbol,
    sequence_expansions,
    metadata,
)
    return (;
        source = :quillen_supplied_local_evidence_patch_assembly,
        local_count = length(candidate.local_certificates),
        raw_denominators = candidate.raw_denominators,
        exponent = solver_result.exponent,
        powered_denominators = solver_result.powered_denominators,
        coverage_multipliers = solver_result.coverage_multipliers,
        coverage_sum = solver_result.coverage_sum,
        substitution_chain_replay_metadata = substitution_chain.replay_metadata,
        base_term_policy = base_term_policy,
        sequence_expansion_metadata = [
            expansion.replay_metadata for expansion in sequence_expansions
        ],
        metadata = metadata,
    )
end

function _same_quillen_sequence_expansions(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        left[idx].local_index == right[idx].local_index || return false
        _same_quillen_denominator_cover_solver_result_data(
            left[idx].solver_result,
            right[idx].solver_result,
        ) || return false
        left[idx].cover_term == right[idx].cover_term || return false
        _same_quillen_factors(
            left[idx].global_elementary_factors,
            right[idx].global_elementary_factors,
        ) || return false
        left[idx].replay_metadata == right[idx].replay_metadata || return false
        _same_quillen_sequence_expansion_verification(
            left[idx].verification,
            right[idx].verification,
        ) || return false
    end
    return true
end

function _same_quillen_supplied_patch_verification(
    left::QuillenSuppliedEvidencePatchAssemblyVerification,
    right::QuillenSuppliedEvidencePatchAssemblyVerification,
)::Bool
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
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _same_quillen_denominator_cover_candidate_data(
    left::QuillenDenominatorCoverCandidate,
    right::QuillenDenominatorCoverCandidate,
)::Bool
    return left.original_input == right.original_input &&
           left.ring == right.ring &&
           left.size == right.size &&
           left.selected_variable == right.selected_variable &&
           left.local_certificates == right.local_certificates &&
           left.raw_denominators == right.raw_denominators &&
           _same_quillen_local_denominator_supports(
               left.local_supports,
               right.local_supports,
           ) &&
           left.replay_metadata == right.replay_metadata &&
           _same_quillen_denominator_cover_candidate_verification(
               left.verification,
               right.verification,
           )
end

function _same_quillen_denominator_cover_candidate_data(left, right)::Bool
    return left === nothing && right === nothing
end

function _same_quillen_denominator_cover_solver_result_data(
    left::QuillenDenominatorCoverSolverResult,
    right::QuillenDenominatorCoverSolverResult,
)::Bool
    return _same_quillen_denominator_cover_candidate_data(
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
end

function _quillen_supplied_base_term_policy(base_term_policy, base_term_factors)
    if base_term_policy === nothing
        base_term_factors === nothing &&
            throw(ArgumentError("Quillen supplied evidence patch assembly requires supplied A(0) factors or base_term_policy = :trivial or :already_handled"))
        return :supplied
    end
    policy = Symbol(base_term_policy)
    policy in (:supplied, :trivial, :already_handled) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly has unsupported base-term policy"))
    policy == :supplied && base_term_factors === nothing &&
        throw(ArgumentError("Quillen supplied evidence patch assembly requires supplied A(0) factors for base_term_policy = :supplied"))
    return policy
end

function _quillen_supplied_base_term_factors(R, n::Int, base_term_factors)
    factor_type = typeof(identity_matrix(R, n))
    base_term_factors === nothing && return factor_type[]
    return [
        _quillen_local_require_factor_matrix(factor, R, n, "base-term factor")
        for factor in collect(base_term_factors)
    ]
end

function _quillen_supplied_base_term_factor_is_elementary(factor, R, n::Int)::Bool
    for row in 1:n
        factor[row, row] == one(R) || return false
    end
    nonzero_offdiagonal_count = 0
    for row in 1:n, col in 1:n
        row == col && continue
        factor[row, col] == zero(R) && continue
        nonzero_offdiagonal_count += 1
        nonzero_offdiagonal_count <= 1 || return false
    end
    return true
end

function _quillen_supplied_base_term_factors_are_elementary(factors, R, n::Int)::Bool
    return all(
        factor -> _quillen_supplied_base_term_factor_is_elementary(factor, R, n),
        factors,
    )
end

function _quillen_supplied_base_term_ok(policy::Symbol, base_term, factors, product, R, n::Int)
    policy == :supplied &&
        return !isempty(factors) &&
               _quillen_supplied_base_term_factors_are_elementary(factors, R, n) &&
               product == base_term
    policy == :trivial && return isempty(factors) && base_term == identity_matrix(R, n)
    policy == :already_handled && return isempty(factors)
    return false
end

function replay_quillen_supplied_evidence_patch(
    patch::QuillenSuppliedEvidencePatchAssembly,
)
    R = _require_quillen_denominator_cover_ring(patch.ring)
    n = patch.size
    _quillen_local_require_factor_matrix(
        patch.original_input,
        R,
        n,
        "supplied evidence original input",
    )
    selected = _require_substitution_generator(R, patch.substitution_variable)
    local_certificates = patch.local_certificates
    local_certificates_ok = !isempty(local_certificates) &&
        all(verify_quillen_local_factor_sequence_certificate, local_certificates)
    expected_candidate = extract_quillen_denominator_cover_candidate(local_certificates)
    denominator_candidate_ok =
        verify_quillen_denominator_cover_candidate(patch.denominator_candidate)
    denominator_candidate_matches =
        denominator_candidate_ok &&
        _same_quillen_denominator_cover_candidate_data(
            patch.denominator_candidate,
            expected_candidate,
        )
    solver_result_ok =
        verify_quillen_denominator_cover_solver_result(patch.solver_result)
    solver_source_candidate_ok =
        solver_result_ok &&
        patch.solver_result.source_candidate isa QuillenDenominatorCoverCandidate &&
        _same_quillen_denominator_cover_candidate_data(
            patch.solver_result.source_candidate,
            patch.denominator_candidate,
        )
    cover_certificate_ok =
        verify_quillen_denominator_cover(patch.cover_certificate) &&
        _same_quillen_cover_certificate_data(
            patch.cover_certificate,
            patch.solver_result.cover_certificate,
        )
    substitution_chain_ok =
        verify_quillen_patch_substitution_chain(patch.substitution_chain)
    substitution_chain_matches =
        substitution_chain_ok &&
        patch.substitution_chain.original_matrix == patch.original_input &&
        patch.substitution_chain.selected_variable == selected &&
        _same_quillen_denominator_cover_solver_result_data(
            patch.substitution_chain.solver_result,
            patch.solver_result,
        )
    base_term_factors = _quillen_supplied_base_term_factors(
        R,
        n,
        patch.base_term_factors,
    )
    base_term_product = _quillen_product(R, n, base_term_factors)
    base_term_ok =
        patch.base_term == patch.substitution_chain.base_term &&
        patch.base_term_product == base_term_product &&
        _quillen_supplied_base_term_ok(
            patch.base_term_policy,
            patch.base_term,
            base_term_factors,
            base_term_product,
            R,
            n,
        )
    sequence_expansions = [
        quillen_sequence_contribution_expansion(certificate, patch.solver_result, index)
        for (index, certificate) in enumerate(local_certificates)
    ]
    sequence_expansions_ok =
        all(verify_quillen_sequence_contribution_expansion, patch.sequence_expansions) &&
        _same_quillen_sequence_expansions(
            patch.sequence_expansions,
            sequence_expansions,
        )
    factor_type = typeof(identity_matrix(R, n))
    sequence_elementary_factors = factor_type[]
    for expansion in sequence_expansions
        append!(sequence_elementary_factors, expansion.global_elementary_factors)
    end
    global_elementary_factors = copy(base_term_factors)
    append!(global_elementary_factors, sequence_elementary_factors)
    global_elementary_factors_ok =
        _same_quillen_factors(
            patch.sequence_elementary_factors,
            sequence_elementary_factors,
        ) &&
        _same_quillen_factors(
            patch.global_elementary_factors,
            global_elementary_factors,
        )
    product = _quillen_product(R, n, global_elementary_factors)
    target = _quillen_local_require_factor_matrix(
        patch.target,
        R,
        n,
        "supplied evidence target",
    )
    product_ok =
        global_elementary_factors_ok &&
        patch.product == product &&
        product == target
    target_ok = target == patch.original_input
    replay_metadata = _quillen_supplied_patch_metadata(
        patch.denominator_candidate,
        patch.solver_result,
        patch.substitution_chain,
        patch.base_term_policy,
        sequence_expansions,
        patch.replay_metadata.metadata,
    )
    replay_metadata_ok = patch.replay_metadata == replay_metadata
    overall_ok =
        local_certificates_ok &&
        denominator_candidate_ok &&
        denominator_candidate_matches &&
        solver_result_ok &&
        solver_source_candidate_ok &&
        cover_certificate_ok &&
        substitution_chain_ok &&
        substitution_chain_matches &&
        base_term_ok &&
        sequence_expansions_ok &&
        global_elementary_factors_ok &&
        product_ok &&
        target_ok &&
        replay_metadata_ok
    return QuillenSuppliedEvidencePatchAssemblyVerification(
        local_certificates_ok,
        denominator_candidate_ok,
        denominator_candidate_matches,
        solver_result_ok,
        solver_source_candidate_ok,
        cover_certificate_ok,
        substitution_chain_ok,
        substitution_chain_matches,
        base_term_ok,
        sequence_expansions_ok,
        global_elementary_factors,
        global_elementary_factors_ok,
        product,
        product_ok,
        target,
        target_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function verify_quillen_patch(patch::QuillenSuppliedEvidencePatchAssembly)::Bool
    try
        replay = replay_quillen_supplied_evidence_patch(patch)
        return replay.overall_ok &&
               _same_quillen_supplied_patch_verification(patch.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function assemble_quillen_patch_from_local_evidence(
    A,
    selected_variable,
    local_certificates;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
    substitution_chain = nothing,
    base_term_policy = nothing,
    base_term_factors = nothing,
    metadata = (;),
)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    n = _require_square_matrix(A, "supplied evidence original input")
    selected = _require_substitution_generator(R, selected_variable)
    certificates = collect(local_certificates)
    isempty(certificates) &&
        throw(ArgumentError("Quillen supplied evidence patch assembly requires at least one local sequence certificate"))
    all(verify_quillen_local_factor_sequence_certificate, certificates) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly requires verified local sequence certificates"))
    all(certificate -> certificate.original_input == A, certificates) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly requires local evidence for the input matrix"))
    candidate = extract_quillen_denominator_cover_candidate(certificates)
    candidate.ring == R &&
        candidate.size == n &&
        candidate.selected_variable == selected ||
        throw(ArgumentError("Quillen supplied evidence patch assembly candidate context does not match the input"))
    solver_result = solve_quillen_denominator_cover(
        candidate;
        max_exponent,
        exponent,
        coverage_multipliers,
        supplied_multipliers,
    )
    chain = substitution_chain === nothing ?
        quillen_patch_substitution_chain(
            A,
            selected,
            solver_result;
            metadata = merge(
                (; source = :quillen_supplied_evidence_patch_assembly),
                metadata,
            ),
        ) :
        substitution_chain
    verify_quillen_patch_substitution_chain(chain) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly substitution chain does not replay"))
    chain.original_matrix == A && chain.selected_variable == selected ||
        throw(ArgumentError("Quillen supplied evidence patch assembly substitution chain context does not match"))
    _same_quillen_denominator_cover_solver_result_data(
        chain.solver_result,
        solver_result,
    ) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly substitution chain solver data does not match"))
    policy = _quillen_supplied_base_term_policy(base_term_policy, base_term_factors)
    normalized_base_factors = _quillen_supplied_base_term_factors(
        R,
        n,
        base_term_factors,
    )
    base_product = _quillen_product(R, n, normalized_base_factors)
    _quillen_supplied_base_term_ok(
        policy,
        chain.base_term,
        normalized_base_factors,
        base_product,
        R,
        n,
    ) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly base-term evidence is missing or does not replay"))
    sequence_expansions = [
        quillen_sequence_contribution_expansion(certificate, solver_result, index)
        for (index, certificate) in enumerate(certificates)
    ]
    factor_type = typeof(identity_matrix(R, n))
    sequence_elementary_factors = factor_type[]
    for expansion in sequence_expansions
        append!(sequence_elementary_factors, expansion.global_elementary_factors)
    end
    global_elementary_factors = copy(normalized_base_factors)
    append!(global_elementary_factors, sequence_elementary_factors)
    product = _quillen_product(R, n, global_elementary_factors)
    product == A ||
        throw(ArgumentError("Quillen supplied evidence patch assembly factors do not multiply to the input matrix"))
    replay_metadata = _quillen_supplied_patch_metadata(
        candidate,
        solver_result,
        chain,
        policy,
        sequence_expansions,
        metadata,
    )
    provisional = QuillenSuppliedEvidencePatchAssembly(
        R,
        n,
        selected,
        A,
        certificates,
        candidate,
        solver_result,
        solver_result.cover_certificate,
        chain,
        policy,
        chain.base_term,
        normalized_base_factors,
        base_product,
        sequence_expansions,
        sequence_elementary_factors,
        global_elementary_factors,
        product,
        A,
        replay_metadata,
        QuillenSuppliedEvidencePatchAssemblyVerification(
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            factor_type[],
            false,
            identity_matrix(R, n),
            false,
            A,
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_quillen_supplied_evidence_patch(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen supplied evidence patch assembly data does not replay"))
    return QuillenSuppliedEvidencePatchAssembly(
        provisional.ring,
        provisional.size,
        provisional.substitution_variable,
        provisional.original_input,
        provisional.local_certificates,
        provisional.denominator_candidate,
        provisional.solver_result,
        provisional.cover_certificate,
        provisional.substitution_chain,
        provisional.base_term_policy,
        provisional.base_term,
        provisional.base_term_factors,
        provisional.base_term_product,
        provisional.sequence_expansions,
        provisional.sequence_elementary_factors,
        provisional.global_elementary_factors,
        provisional.product,
        provisional.target,
        provisional.replay_metadata,
        verification,
    )
end

struct QuillenMurthyAdapterConsumptionVerification
    adapters_ok::Bool
    context_ok::Bool
    local_sequences_ok::Bool
    patch_ok::Bool
    patch_alignment_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenMurthyAdapterConsumption
    original_input
    ring
    size::Int
    selected_variable
    murthy_adapters::Vector{MurthyQuillenLocalAdapter}
    local_sequence_certificates::Vector{QuillenLocalFactorSequenceCertificate}
    patch::QuillenSuppliedEvidencePatchAssembly
    replay_metadata
    verification::QuillenMurthyAdapterConsumptionVerification
end

function _same_quillen_murthy_adapter_consumption_verification(
    left::QuillenMurthyAdapterConsumptionVerification,
    right::QuillenMurthyAdapterConsumptionVerification,
)::Bool
    return left.adapters_ok == right.adapters_ok &&
           left.context_ok == right.context_ok &&
           left.local_sequences_ok == right.local_sequences_ok &&
           left.patch_ok == right.patch_ok &&
           left.patch_alignment_ok == right.patch_alignment_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _quillen_murthy_adapter_vector(adapters)
    collected = collect(adapters)
    isempty(collected) &&
        throw(ArgumentError("Murthy adapter consumption requires at least one adapter"))
    normalized = MurthyQuillenLocalAdapter[]
    for adapter in collected
        adapter isa MurthyQuillenLocalAdapter ||
            throw(ArgumentError("Murthy adapter consumption requires MurthyQuillenLocalAdapter records"))
        push!(normalized, adapter)
    end
    return normalized
end

function _same_sl3_local_factor_replay(
    left::SL3LocalElementaryFactorReplay,
    right::SL3LocalElementaryFactorReplay,
)::Bool
    return left.target == right.target &&
           left.factors == right.factors &&
           left.selected_variable == right.selected_variable &&
           left.mode == right.mode &&
           left.denominator_product == right.denominator_product &&
           left.cleared_product == right.cleared_product &&
           left.materialized_factors == right.materialized_factors
end

function _quillen_murthy_adapter_patch_metadata(adapters, local_sequences, metadata)
    return (;
        source = :quillen_murthy_adapter_consumption,
        adapter_count = length(adapters),
        murthy_adapter_metadata = [adapter.replay_metadata for adapter in adapters],
        local_sequence_metadata = [
            certificate.replay_metadata for certificate in local_sequences
        ],
        metadata = metadata,
    )
end

function _quillen_murthy_adapter_consumption_metadata(
    A,
    selected_variable,
    adapters,
    local_sequences,
    patch::QuillenSuppliedEvidencePatchAssembly,
    metadata,
)
    return (;
        source = :quillen_murthy_adapter_consumption,
        adapter_count = length(adapters),
        selected_variable = selected_variable,
        original_input = A,
        murthy_adapter_metadata = [adapter.replay_metadata for adapter in adapters],
        local_sequence_metadata = [
            certificate.replay_metadata for certificate in local_sequences
        ],
        patch_metadata = patch.replay_metadata,
        metadata = metadata,
    )
end

function _quillen_murthy_adapter_context_check(
    adapter::MurthyQuillenLocalAdapter,
    A,
    R,
    n::Int,
    selected_variable,
)
    _verify_murthy_quillen_local_adapter(adapter) ||
        throw(ArgumentError("Murthy adapter consumption requires verified Murthy Quillen adapters"))
    adapter.ring == R ||
        throw(ArgumentError("Murthy adapter consumption adapter ring does not match the input"))
    adapter.size == n ||
        throw(ArgumentError("Murthy adapter consumption adapter size does not match the input"))
    adapter.original_input == A ||
        throw(ArgumentError("Murthy adapter consumption adapter original input does not match the input"))
    adapter.selected_variable == selected_variable ||
        throw(ArgumentError("Murthy adapter consumption selected variable mismatch"))
    adapter.murthy_certificate.target == A ||
        throw(ArgumentError("Murthy adapter consumption Murthy target does not match the input"))
    adapter.murthy_certificate.selected_variable == selected_variable ||
        throw(ArgumentError("Murthy adapter consumption Murthy selected variable mismatch"))

    replay = adapter.local_factor_replay
    verify_sl3_local_elementary_factor_replay(replay) ||
        throw(ArgumentError("Murthy adapter consumption local factor replay does not verify"))
    expected_replay = _murthy_quillen_local_replay(adapter.murthy_certificate)
    _same_sl3_local_factor_replay(replay, expected_replay) ||
        throw(ArgumentError("Murthy adapter consumption local factor replay does not match the Murthy certificate"))
    replay.target == A ||
        throw(ArgumentError("Murthy adapter consumption local factor replay target mismatch"))
    replay.selected_variable == selected_variable ||
        throw(ArgumentError("Murthy adapter consumption local factor replay selected variable mismatch"))
    replay.denominator_product == prod(
        (factor.denominator for factor in replay.factors);
        init = one(R),
    ) ||
        throw(ArgumentError("Murthy adapter consumption denominator product mismatch"))
    adapter.materialized_factors == replay.materialized_factors ||
        throw(ArgumentError("Murthy adapter consumption materialized factor mismatch"))
    expected_metadata = _murthy_quillen_local_replay_metadata(
        adapter.murthy_certificate,
        replay,
        adapter.mode,
        adapter.witness_metadata,
    )
    adapter.replay_metadata == expected_metadata ||
        throw(ArgumentError("Murthy adapter consumption replay metadata mismatch"))
    return true
end

function _quillen_murthy_sequence_has_denominator_provenance(
    certificate::QuillenLocalFactorSequenceCertificate,
)::Bool
    expected_denominator_data = _quillen_denominator_data(
        certificate.normalized_local_contributions,
    )
    return certificate.raw_denominators == [
               factor.denominator for factor in certificate.factors
           ] &&
           certificate.product_denominator ==
               prod(certificate.raw_denominators; init = one(certificate.ring)) &&
           _same_quillen_denominator_data(
               certificate.verification.denominator_data,
               expected_denominator_data,
           ) &&
           hasproperty(certificate.replay_metadata, :denominator_data) &&
           _same_quillen_denominator_data(
               certificate.replay_metadata.denominator_data,
               expected_denominator_data,
           ) &&
           hasproperty(certificate.replay_metadata, :factor_provenance)
end

function _quillen_murthy_factor_provenance_ok(
    certificate::QuillenLocalFactorSequenceCertificate,
)::Bool
    return all(certificate.factors) do factor
        present, source = _quillen_local_sequence_provenance_field(
            factor.provenance,
            :source,
        )
        present && source == :murthy_local_sl3
    end
end

function _same_quillen_local_factor_sequence_certificate_data(
    left::QuillenLocalFactorSequenceCertificate,
    right::QuillenLocalFactorSequenceCertificate,
)::Bool
    return left.original_input == right.original_input &&
           left.ring == right.ring &&
           left.size == right.size &&
           left.selected_variable == right.selected_variable &&
           _same_quillen_local_elementary_factors(left.factors, right.factors) &&
           left.raw_denominators == right.raw_denominators &&
           left.product_denominator == right.product_denominator &&
           left.local_product == right.local_product &&
           left.local_correction == right.local_correction &&
           _same_quillen_local_contributions(
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

function _same_quillen_local_factor_sequence_certificates(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        _same_quillen_local_factor_sequence_certificate_data(
            left[idx],
            right[idx],
        ) || return false
    end
    return true
end

function _quillen_murthy_rebuilt_local_factor_sequence(
    adapter::MurthyQuillenLocalAdapter,
)
    return _murthy_quillen_local_factor_sequence(
        adapter.original_input,
        adapter.selected_variable,
        adapter.local_factor_replay,
        adapter.witness_metadata,
        (; source = :murthy_quillen_local_adapter),
    )
end

function _quillen_murthy_adapter_sequences(A, R, n::Int, selected_variable, adapters)
    sequences = QuillenLocalFactorSequenceCertificate[]
    for adapter in adapters
        _quillen_murthy_adapter_context_check(adapter, A, R, n, selected_variable)
        cached_sequence = _murthy_quillen_local_factor_sequence_certificate(adapter)
        expected_sequence = _quillen_murthy_rebuilt_local_factor_sequence(adapter)
        _same_quillen_local_factor_sequence_certificate_data(
            cached_sequence,
            expected_sequence,
        ) ||
            throw(ArgumentError("Murthy adapter consumption cached local sequence does not match the Murthy replay"))
        sequence = expected_sequence
        verify_quillen_local_factor_sequence_certificate(sequence) ||
            throw(ArgumentError("Murthy adapter consumption converted local sequence does not replay"))
        sequence.original_input == A ||
            throw(ArgumentError("Murthy adapter consumption converted sequence original input mismatch"))
        sequence.ring == R && sequence.size == n ||
            throw(ArgumentError("Murthy adapter consumption converted sequence context mismatch"))
        sequence.selected_variable == selected_variable ||
            throw(ArgumentError("Murthy adapter consumption converted sequence selected variable mismatch"))
        sequence.local_product == adapter.local_product &&
            sequence.local_correction == adapter.local_correction ||
            throw(ArgumentError("Murthy adapter consumption converted sequence product mismatch"))
        _quillen_murthy_sequence_has_denominator_provenance(sequence) ||
            throw(ArgumentError("Murthy adapter consumption converted sequence lacks denominator provenance"))
        _quillen_murthy_factor_provenance_ok(sequence) ||
            throw(ArgumentError("Murthy adapter consumption converted sequence lacks Murthy factor provenance"))
        push!(sequences, sequence)
    end
    return sequences
end

function quillen_local_sequences_from_murthy_adapters(A, selected_variable, adapters)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    n = _require_square_matrix(A, "Murthy adapter consumption original input")
    n == 3 || throw(ArgumentError("Murthy adapter consumption requires a 3x3 input"))
    selected = _require_substitution_generator(R, selected_variable)
    collected = _quillen_murthy_adapter_vector(adapters)
    return _quillen_murthy_adapter_sequences(A, R, n, selected, collected)
end

function assemble_quillen_patch_from_murthy_adapters(
    A,
    selected_variable,
    adapters;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
    substitution_chain = nothing,
    base_term_policy = nothing,
    base_term_factors = nothing,
    metadata = (;),
)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    n = _require_square_matrix(A, "Murthy adapter consumption original input")
    n == 3 || throw(ArgumentError("Murthy adapter consumption requires a 3x3 input"))
    selected = _require_substitution_generator(R, selected_variable)
    collected = _quillen_murthy_adapter_vector(adapters)
    sequences = _quillen_murthy_adapter_sequences(A, R, n, selected, collected)
    patch_metadata = _quillen_murthy_adapter_patch_metadata(
        collected,
        sequences,
        metadata,
    )
    return assemble_quillen_patch_from_local_evidence(
        A,
        selected,
        sequences;
        max_exponent,
        exponent,
        coverage_multipliers,
        supplied_multipliers,
        substitution_chain,
        base_term_policy,
        base_term_factors,
        metadata = patch_metadata,
    )
end

function _quillen_murthy_consumption_verification(
    A,
    R,
    n::Int,
    selected_variable,
    adapters,
    local_sequences,
    patch::QuillenSuppliedEvidencePatchAssembly,
    replay_metadata,
)
    adapters_ok = false
    context_ok = false
    local_sequences_ok = false
    expected_sequences = QuillenLocalFactorSequenceCertificate[]
    try
        adapters_ok = all(adapter -> _verify_murthy_quillen_local_adapter(adapter), adapters)
        context_ok = adapters_ok && all(adapters) do adapter
            adapter.ring == R &&
                adapter.size == n &&
                adapter.original_input == A &&
                adapter.selected_variable == selected_variable
        end
        expected_sequences = context_ok ?
            _quillen_murthy_adapter_sequences(A, R, n, selected_variable, adapters) :
            QuillenLocalFactorSequenceCertificate[]
        local_sequences_ok =
            context_ok &&
            _same_quillen_local_factor_sequence_certificates(
                local_sequences,
                expected_sequences,
            ) &&
            all(verify_quillen_local_factor_sequence_certificate, local_sequences)
    catch err
        err isa InterruptException && rethrow()
        local_sequences_ok = false
    end
    patch_ok = verify_quillen_patch(patch)
    metadata = hasproperty(replay_metadata, :metadata) ?
        replay_metadata.metadata :
        nothing
    expected_patch_metadata = _quillen_murthy_adapter_patch_metadata(
        adapters,
        local_sequences,
        metadata,
    )
    patch_metadata_ok =
        patch_ok &&
        hasproperty(patch.replay_metadata, :metadata) &&
        patch.replay_metadata.metadata == expected_patch_metadata
    patch_alignment_ok =
        patch_ok &&
        local_sequences_ok &&
        patch.original_input == A &&
        patch.ring == R &&
        patch.size == n &&
        patch.substitution_variable == selected_variable &&
        _same_quillen_local_factor_sequence_certificates(
            patch.local_certificates,
            local_sequences,
        ) &&
        patch_metadata_ok
    expected_metadata = _quillen_murthy_adapter_consumption_metadata(
        A,
        selected_variable,
        adapters,
        local_sequences,
        patch,
        metadata,
    )
    replay_metadata_ok = replay_metadata == expected_metadata
    overall_ok =
        adapters_ok &&
        context_ok &&
        local_sequences_ok &&
        patch_ok &&
        patch_alignment_ok &&
        replay_metadata_ok
    return QuillenMurthyAdapterConsumptionVerification(
        adapters_ok,
        context_ok,
        local_sequences_ok,
        patch_ok,
        patch_alignment_ok,
        expected_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function replay_quillen_murthy_adapter_consumption(
    consumption::QuillenMurthyAdapterConsumption,
)
    R = _require_quillen_denominator_cover_ring(consumption.ring)
    _quillen_local_require_factor_matrix(
        consumption.original_input,
        R,
        consumption.size,
        "Murthy adapter consumption original input",
    )
    selected = _require_substitution_generator(R, consumption.selected_variable)
    return _quillen_murthy_consumption_verification(
        consumption.original_input,
        R,
        consumption.size,
        selected,
        consumption.murthy_adapters,
        consumption.local_sequence_certificates,
        consumption.patch,
        consumption.replay_metadata,
    )
end

function verify_quillen_murthy_adapter_consumption(consumption)::Bool
    try
        replay = replay_quillen_murthy_adapter_consumption(consumption)
        return replay.overall_ok &&
               _same_quillen_murthy_adapter_consumption_verification(
                   consumption.verification,
                   replay,
               )
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function consume_murthy_quillen_adapters_for_patch(
    A,
    selected_variable,
    adapters;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
    substitution_chain = nothing,
    base_term_policy = nothing,
    base_term_factors = nothing,
    metadata = (;),
)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    n = _require_square_matrix(A, "Murthy adapter consumption original input")
    n == 3 || throw(ArgumentError("Murthy adapter consumption requires a 3x3 input"))
    selected = _require_substitution_generator(R, selected_variable)
    collected = _quillen_murthy_adapter_vector(adapters)
    sequences = _quillen_murthy_adapter_sequences(A, R, n, selected, collected)
    patch = assemble_quillen_patch_from_local_evidence(
        A,
        selected,
        sequences;
        max_exponent,
        exponent,
        coverage_multipliers,
        supplied_multipliers,
        substitution_chain,
        base_term_policy,
        base_term_factors,
        metadata = _quillen_murthy_adapter_patch_metadata(
            collected,
            sequences,
            metadata,
        ),
    )
    replay_metadata = _quillen_murthy_adapter_consumption_metadata(
        A,
        selected,
        collected,
        sequences,
        patch,
        metadata,
    )
    placeholder = QuillenMurthyAdapterConsumptionVerification(
        false,
        false,
        false,
        false,
        false,
        replay_metadata,
        false,
        false,
    )
    provisional = QuillenMurthyAdapterConsumption(
        A,
        R,
        n,
        selected,
        collected,
        sequences,
        patch,
        replay_metadata,
        placeholder,
    )
    verification = replay_quillen_murthy_adapter_consumption(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Murthy adapter consumption data does not replay"))
    return QuillenMurthyAdapterConsumption(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.murthy_adapters,
        provisional.local_sequence_certificates,
        provisional.patch,
        provisional.replay_metadata,
        verification,
    )
end

struct QuillenGlobalPatchAssemblyVerification
    cover_certificate_ok::Bool
    local_certificates_ok::Bool
    normalized_contributions_ok::Bool
    local_count_ok::Bool
    local_alignment_ok::Bool
    cover_alignment_ok::Bool
    normalized_input_ok::Bool
    selected_variable_ok::Bool
    denominator_data_ok::Bool
    coverage_sum
    coverage_ok::Bool
    global_elementary_factors::Vector
    global_elementary_factors_ok::Bool
    product
    product_ok::Bool
    target
    target_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenGlobalPatchAssembly
    ring
    size::Int
    substitution_variable
    original_input
    cover_certificate::QuillenDenominatorCoverCertificate
    denominator_data::Vector{QuillenDenominatorData}
    local_certificates::Vector{QuillenLocalRealizationCertificate}
    normalized_local_contributions::Vector{QuillenLocalContributionNormalization}
    global_elementary_factors::Vector
    patched_product
    target
    replay_metadata
    verification::QuillenGlobalPatchAssemblyVerification
end

function _same_quillen_local_certificate_data(
    left::QuillenLocalRealizationCertificate,
    right::QuillenLocalRealizationCertificate,
)::Bool
    return left.original_input == right.original_input &&
           left.ring == right.ring &&
           left.size == right.size &&
           left.selected_variable == right.selected_variable &&
           left.local_certificate.indices == right.local_certificate.indices &&
           left.local_certificate.denominators == right.local_certificate.denominators &&
           left.denominator == right.denominator &&
           left.coverage_multiplier == right.coverage_multiplier &&
           left.correction == right.correction &&
           _same_quillen_factors(left.factors, right.factors) &&
           left.local_product == right.local_product &&
           left.local_correction == right.local_correction &&
           left.patched_substitution_witness == right.patched_substitution_witness &&
           left.witness_metadata == right.witness_metadata &&
           left.verification == right.verification
end

function _same_quillen_cover_certificate_data(
    left::QuillenDenominatorCoverCertificate,
    right::QuillenDenominatorCoverCertificate,
)::Bool
    return left.ring == right.ring &&
           left.denominators == right.denominators &&
           left.coverage_multipliers == right.coverage_multipliers &&
           left.coverage_sum == right.coverage_sum &&
           _same_quillen_denominator_cover_verification(left.verification, right.verification)
end

function _quillen_global_patch_replay_metadata(
    cover::QuillenDenominatorCoverCertificate,
    local_certificates,
    normalized_contributions,
)
    return (;
        fixture_count = length(local_certificates),
        normalized_count = length(normalized_contributions),
        cover_denominator_count = length(cover.denominators),
        cover_coverage_sum = cover.coverage_sum,
        local_witness_metadata = [certificate.witness_metadata for certificate in local_certificates],
        normalized_cover_indices = [normalized.cover_index for normalized in normalized_contributions],
    )
end

function _same_quillen_global_patch_verification(
    left::QuillenGlobalPatchAssemblyVerification,
    right::QuillenGlobalPatchAssemblyVerification,
)::Bool
    return left.cover_certificate_ok == right.cover_certificate_ok &&
           left.local_certificates_ok == right.local_certificates_ok &&
           left.normalized_contributions_ok == right.normalized_contributions_ok &&
           left.local_count_ok == right.local_count_ok &&
           left.local_alignment_ok == right.local_alignment_ok &&
           left.cover_alignment_ok == right.cover_alignment_ok &&
           left.normalized_input_ok == right.normalized_input_ok &&
           left.selected_variable_ok == right.selected_variable_ok &&
           left.denominator_data_ok == right.denominator_data_ok &&
           left.coverage_sum == right.coverage_sum &&
           left.coverage_ok == right.coverage_ok &&
           _same_quillen_factors(
               left.global_elementary_factors,
               right.global_elementary_factors,
           ) &&
           left.global_elementary_factors_ok == right.global_elementary_factors_ok &&
           left.product == right.product &&
           left.product_ok == right.product_ok &&
           left.target == right.target &&
           left.target_ok == right.target_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _quillen_global_target_matrix(target; ring, size::Int, label::AbstractString)
    if target isa QuillenElementaryCorrection
        correction = _normalize_quillen_contribution(
            QuillenLocalContribution(
                LocalCertificate([target.row, target.col], [one(ring), one(ring)]),
                one(ring),
                one(ring),
                target,
            ),
            ring,
            size,
        ).correction
        return elementary_matrix(size, correction.row, correction.col, correction.entry, ring)
    end
    _quillen_local_require_factor_matrix(target, ring, size, label)
    return target
end

function replay_deterministic_quillen_patch(patch::QuillenGlobalPatchAssembly)
    R = _require_supported_quillen_ring(patch.ring)
    n = patch.size
    n >= 2 || throw(ArgumentError("patch size must be at least 2"))
    selected = _require_substitution_generator(R, patch.substitution_variable)
    target = _quillen_global_target_matrix(
        patch.target;
        ring = R,
        size = n,
        label = "target",
    )

    cover = patch.cover_certificate
    local_certificates = patch.local_certificates
    normalized = patch.normalized_local_contributions
    cover_certificate_ok = verify_quillen_denominator_cover(cover)
    local_certificates_ok = all(verify_quillen_local_certificate, local_certificates)
    normalized_contributions_ok =
        all(verify_quillen_local_contribution_normalization, normalized)
    local_count_ok = length(local_certificates) == length(normalized)

    local_alignment_ok = local_count_ok && all(eachindex(local_certificates)) do idx
        _same_quillen_local_certificate_data(
            normalized[idx].local_certificate,
            local_certificates[idx],
        )
    end
    cover_alignment_ok =
        _same_quillen_cover_certificate_data(cover, patch.cover_certificate) &&
        all(normalized) do item
            _same_quillen_cover_certificate_data(item.cover_certificate, cover)
        end
    normalized_input_ok = all(normalized) do item
        item.original_input == patch.original_input &&
            item.local_certificate.original_input == patch.original_input
    end
    selected_variable_ok =
        all(normalized) do item
            item.selected_variable == selected &&
                item.local_certificate.selected_variable == selected
        end

    expected_denominator_data = _quillen_denominator_data([
        item.local_contribution for item in normalized
    ])
    denominator_data_ok =
        _same_quillen_denominator_data(patch.denominator_data, expected_denominator_data)
    coverage_sum = _quillen_coverage_sum(R, patch.denominator_data)
    coverage_ok =
        denominator_data_ok && cover_certificate_ok && coverage_sum == one(R) &&
        coverage_sum == cover.coverage_sum
    global_elementary_factors =
        reduce(vcat, [item.weighted_global_elementary_factors for item in normalized]; init = Any[])
    global_elementary_factors_ok = _same_quillen_factors(
        patch.global_elementary_factors,
        global_elementary_factors,
    )
    product = _quillen_product(R, n, patch.global_elementary_factors)
    product_ok = global_elementary_factors_ok && patch.patched_product == product
    target_ok = product == target
    replay_metadata =
        _quillen_global_patch_replay_metadata(cover, local_certificates, normalized)
    replay_metadata_ok = patch.replay_metadata == replay_metadata
    overall_ok =
        cover_certificate_ok &&
        local_certificates_ok &&
        normalized_contributions_ok &&
        local_count_ok &&
        local_alignment_ok &&
        cover_alignment_ok &&
        normalized_input_ok &&
        selected_variable_ok &&
        denominator_data_ok &&
        coverage_ok &&
        global_elementary_factors_ok &&
        product_ok &&
        target_ok &&
        replay_metadata_ok

    return QuillenGlobalPatchAssemblyVerification(
        cover_certificate_ok,
        local_certificates_ok,
        normalized_contributions_ok,
        local_count_ok,
        local_alignment_ok,
        cover_alignment_ok,
        normalized_input_ok,
        selected_variable_ok,
        denominator_data_ok,
        coverage_sum,
        coverage_ok,
        global_elementary_factors,
        global_elementary_factors_ok,
        product,
        product_ok,
        target,
        target_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function verify_quillen_patch(patch::QuillenGlobalPatchAssembly)::Bool
    try
        replay = replay_deterministic_quillen_patch(patch)
        return replay.overall_ok &&
               _same_quillen_global_patch_verification(patch.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function assemble_deterministic_quillen_patch(
    original_input,
    selected_variable,
    local_certificates,
    normalized_local_contributions,
    cover::QuillenDenominatorCoverCertificate;
    target = original_input,
    ring = nothing,
    size = nothing,
)
    certificates = collect(local_certificates)
    normalized = collect(normalized_local_contributions)
    isempty(certificates) &&
        throw(ArgumentError("deterministic Quillen patch assembly requires at least one local certificate"))
    length(certificates) == length(normalized) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires matching local certificate and normalized contribution counts"))
    verify_quillen_denominator_cover(cover) ||
        throw(ArgumentError("Quillen denominator cover certificate does not replay"))
    all(verify_quillen_local_certificate, certificates) ||
        throw(ArgumentError("Quillen local realization certificate does not replay"))
    all(verify_quillen_local_contribution_normalization, normalized) ||
        throw(ArgumentError("Quillen local contribution normalization does not replay"))

    R, n = _quillen_local_input_ring_size(original_input; ring, size)
    cover.ring == R ||
        throw(ArgumentError("deterministic Quillen patch assembly requires cover and target rings to match"))
    selected = _require_substitution_generator(R, selected_variable)
    target_matrix = _quillen_global_target_matrix(target; ring = R, size = n, label = "target")

    all(certificate -> certificate.ring == R && certificate.size == n, certificates) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires local certificate rings and sizes to match the target"))
    all(certificate -> certificate.selected_variable == selected, certificates) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires local certificate variables to match"))
    all(certificate -> certificate.original_input == original_input, certificates) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires local certificate original inputs to match"))

    for idx in eachindex(certificates)
        _same_quillen_local_certificate_data(normalized[idx].local_certificate, certificates[idx]) ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution local certificates to match input order"))
        _same_quillen_cover_certificate_data(normalized[idx].cover_certificate, cover) ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution cover certificates to match"))
        normalized[idx].original_input == original_input ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution original inputs to match"))
        normalized[idx].selected_variable == selected ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution selected variables to match"))
    end

    denominator_data = _quillen_denominator_data([
        item.local_contribution for item in normalized
    ])
    coverage_sum = _quillen_coverage_sum(R, denominator_data)
    coverage_sum == one(R) ||
        throw(ArgumentError("deterministic Quillen patch assembly denominator coverage must sum to one"))

    factor_type = typeof(identity_matrix(R, n))
    global_elementary_factors = factor_type[]
    for item in normalized
        append!(global_elementary_factors, item.weighted_global_elementary_factors)
    end
    patched_product = _quillen_product(R, n, global_elementary_factors)
    patched_product == target_matrix ||
        throw(ArgumentError("deterministic Quillen patch assembly product does not equal the target"))

    replay_metadata =
        _quillen_global_patch_replay_metadata(cover, certificates, normalized)
    provisional = QuillenGlobalPatchAssembly(
        R,
        n,
        selected,
        original_input,
        cover,
        denominator_data,
        certificates,
        normalized,
        global_elementary_factors,
        patched_product,
        target_matrix,
        replay_metadata,
        QuillenGlobalPatchAssemblyVerification(
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            zero(R),
            false,
            Any[],
            false,
            identity_matrix(R, n),
            false,
            target_matrix,
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_deterministic_quillen_patch(provisional)
    verification.overall_ok ||
        throw(ArgumentError("deterministic Quillen patch assembly data does not replay"))
    return QuillenGlobalPatchAssembly(
        provisional.ring,
        provisional.size,
        provisional.substitution_variable,
        provisional.original_input,
        provisional.cover_certificate,
        provisional.denominator_data,
        provisional.local_certificates,
        provisional.normalized_local_contributions,
        provisional.global_elementary_factors,
        provisional.patched_product,
        provisional.target,
        provisional.replay_metadata,
        verification,
    )
end
