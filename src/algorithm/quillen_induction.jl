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
