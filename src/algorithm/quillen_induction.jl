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
