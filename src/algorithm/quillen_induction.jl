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
