# Backend-sensitive scaffolding for the later Quillen patching layer.
struct LocalCertificate
    indices::Vector{Int}
    denominators::Vector

    function LocalCertificate(indices::AbstractVector{<:Integer}, denominators::AbstractVector)
        Base.require_one_based_indexing(indices)
        Base.require_one_based_indexing(denominators)
        length(indices) == length(denominators) ||
            throw(ArgumentError("indices and denominators must have the same length"))
        return new(Int.(collect(indices)), collect(denominators))
    end
end

function common_denominator_factor(entries::AbstractVector)
    Base.require_one_based_indexing(entries)
    isempty(entries) && throw(ArgumentError("entries must not be empty"))

    # Task-7 scaffolding keeps the contract simple: return an exact product that clears
    # all toy denominators, not a normalized least common multiple.
    factor = _denominator_factor(entries[1])
    for idx in 2:length(entries)
        factor *= _denominator_factor(entries[idx])
    end
    return factor
end

function _denominator_factor(entry)
    applicable(denominator, entry) && return denominator(entry)
    return one(parent(entry))
end

struct QuillenDenominatorData
    denominator
    coverage_multiplier
end

struct QuillenElementaryCorrection
    row::Int
    col::Int
    entry
end

struct QuillenLocalContribution
    certificate::LocalCertificate
    denominator
    coverage_multiplier
    correction::QuillenElementaryCorrection
end

struct QuillenPatchVerification
    denominator_data_ok::Bool
    coverage_sum
    coverage_ok::Bool
    product
    target
    product_ok::Bool
end

struct QuillenPatch
    ring
    size::Int
    substitution_variable
    denominator_data::Vector{QuillenDenominatorData}
    local_contributions::Vector{QuillenLocalContribution}
    factors::Vector
    product
    target
    verification::QuillenPatchVerification
end

function _require_supported_quillen_ring(R)
    _is_laurent_polynomial_ring(R) && return R
    try
        gens(R)
        coefficient_ring(R)
        return R
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        throw(ArgumentError("target base ring must be a supported exact polynomial or Laurent polynomial ring"))
    end
end

function _require_quillen_target(target, n::Int)
    n >= 2 || throw(ArgumentError("patch size must be at least 2"))
    size = _require_square_matrix(target, "target correction")
    size == n || throw(DimensionMismatch("target correction size must match the requested patch size"))
    return _require_supported_quillen_ring(base_ring(target))
end

function _require_substitution_generator(R, X)
    coerced = _coerce_into_ring(R, X, "substitution variable")
    any(gen -> gen == coerced, collect(gens(R))) ||
        throw(ArgumentError("substitution variable must be a generator of the target ring"))
    return coerced
end

function _require_elementary_indices(n::Int, row::Int, col::Int)
    1 <= row <= n || throw(ArgumentError("elementary correction row must be between 1 and the patch size"))
    1 <= col <= n || throw(ArgumentError("elementary correction column must be between 1 and the patch size"))
    row != col || throw(ArgumentError("elementary correction row and column must differ"))
    return row, col
end

function _coerced_certificate_denominators(certificate::LocalCertificate, R)
    return [_coerce_into_ring(R, denominator, "certificate denominator") for denominator in certificate.denominators]
end

function _normalize_quillen_contribution(contribution::QuillenLocalContribution, R, n::Int)
    correction = contribution.correction
    row, col = _require_elementary_indices(n, correction.row, correction.col)
    row in contribution.certificate.indices && col in contribution.certificate.indices ||
        throw(ArgumentError("local certificate indices must include the correction row and column"))

    denominator = _coerce_into_ring(R, contribution.denominator, "denominator")
    certificate_denominators = _coerced_certificate_denominators(contribution.certificate, R)
    any(certificate_denominator -> certificate_denominator == denominator, certificate_denominators) ||
        throw(ArgumentError("local contribution denominator must be present in its certificate denominators"))

    coverage_multiplier = _coerce_into_ring(R, contribution.coverage_multiplier, "coverage multiplier")
    entry = _coerce_into_ring(R, correction.entry, "correction entry")
    normalized_correction = QuillenElementaryCorrection(row, col, entry)
    normalized_certificate = LocalCertificate(contribution.certificate.indices, certificate_denominators)
    return QuillenLocalContribution(normalized_certificate, denominator, coverage_multiplier, normalized_correction)
end

function _quillen_denominator_data(local_contributions)
    return [
        QuillenDenominatorData(contribution.denominator, contribution.coverage_multiplier)
        for contribution in local_contributions
    ]
end

function _quillen_coverage_sum(R, denominator_data)
    total = zero(R)
    for data in denominator_data
        total += data.coverage_multiplier * data.denominator
    end
    return total
end

function _quillen_factors(R, n::Int, local_contributions)
    factor_type = typeof(identity_matrix(R, n))
    factors = factor_type[]
    for contribution in local_contributions
        correction = contribution.correction
        weighted_entry = contribution.coverage_multiplier * contribution.denominator * correction.entry
        push!(factors, elementary_matrix(n, correction.row, correction.col, weighted_entry, R))
    end
    return factors
end

function _quillen_product(R, n::Int, factors)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _same_quillen_denominator_data(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        left[idx] == right[idx] || return false
    end
    return true
end

function _quillen_patch_verification(R, n::Int, denominator_data, local_contributions, factors, product, target)
    expected_denominator_data = _quillen_denominator_data(local_contributions)
    denominator_data_ok = _same_quillen_denominator_data(denominator_data, expected_denominator_data)
    coverage_sum = _quillen_coverage_sum(R, denominator_data)
    coverage_ok = coverage_sum == one(R)
    actual_product = _quillen_product(R, n, factors)
    product_ok = actual_product == target && product == actual_product
    return QuillenPatchVerification(denominator_data_ok, coverage_sum, coverage_ok, actual_product, target, product_ok)
end

function construct_quillen_patch(n::Int, X, contributions; target)
    collected = collect(contributions)
    isempty(collected) && throw(ArgumentError("local contributions must be nonempty"))

    R = _require_quillen_target(target, n)
    substitution_variable = _require_substitution_generator(R, X)
    local_contributions = [
        _normalize_quillen_contribution(contribution, R, n)
        for contribution in collected
    ]
    denominator_data = _quillen_denominator_data(local_contributions)
    coverage_sum = _quillen_coverage_sum(R, denominator_data)
    coverage_sum == one(R) || throw(ArgumentError("denominator coverage must sum to one"))

    factors = _quillen_factors(R, n, local_contributions)
    product = _quillen_product(R, n, factors)
    product == target || throw(ArgumentError("constructed Quillen patch does not multiply to the target correction"))
    verification = _quillen_patch_verification(R, n, denominator_data, local_contributions, factors, product, target)
    return QuillenPatch(R, n, substitution_variable, denominator_data, local_contributions, factors, product, target, verification)
end

function verify_quillen_patch(patch::QuillenPatch)::Bool
    try
        verification = _quillen_patch_verification(
            patch.ring,
            patch.size,
            patch.denominator_data,
            patch.local_contributions,
            patch.factors,
            patch.product,
            patch.target,
        )
        return verification.denominator_data_ok && verification.coverage_ok && verification.product_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
