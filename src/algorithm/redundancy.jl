struct SteinbergOptimizationCertificate
    original_factors::Vector
    optimized_factors::Vector
    applied_rewrites::Vector
    comparison_summary
    original_product
    optimized_product
    verification
end

function _require_steinberg_ordinary_polynomial_ring(R)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("Steinberg optimization certificates require an ordinary polynomial ring"))
    R isa MPolyRing ||
        throw(ArgumentError("Steinberg optimization certificates require an ordinary polynomial ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("Steinberg optimization certificates require an exact coefficient ring"))
    coefficient_ring(R) isa Field ||
        throw(ArgumentError("Steinberg optimization certificates require a field-backed polynomial ring"))
    return R
end

function _steinberg_sequence_context(factors, label::AbstractString, expected = nothing)
    collected = collect(factors)
    if isempty(collected)
        expected === nothing &&
            throw(ArgumentError("$(label) factor sequence must be nonempty"))
        return (; factors = collected, n = expected.n, ring = expected.ring, records = Any[])
    end

    if expected === nothing
        n = _require_square_matrix(first(collected), "$(label) factor")
        R = _require_steinberg_ordinary_polynomial_ring(base_ring(first(collected)))
    else
        n = expected.n
        R = expected.ring
    end

    records = Any[]
    for (index, factor) in enumerate(collected)
        factor_n = _require_square_matrix(factor, "$(label) factor[$index]")
        factor_n == n ||
            throw(ArgumentError("$(label) factor[$index] must have the common matrix size"))
        _same_base_ring(base_ring(factor), R) ||
            throw(ArgumentError("$(label) factor[$index] must have the common base ring"))
        record = _canonical_elementary_factor_record(factor)
        _same_base_ring(record.ring, R) ||
            throw(ArgumentError("$(label) factor[$index] canonical record ring mismatch"))
        push!(records, record)
    end

    return (; factors = collected, n, ring = R, records)
end

function _steinberg_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for (index, factor) in enumerate(factors)
        nrows(factor) == n && ncols(factor) == n ||
            throw(ArgumentError("Steinberg factor[$index] must have the common matrix size"))
        _same_base_ring(base_ring(factor), R) ||
            throw(ArgumentError("Steinberg factor[$index] must have the common base ring"))
        product *= factor
    end
    return product
end

function _steinberg_metric_summary(factors)
    return (;
        max_elementary_factor_monomial_degree = max_elementary_factor_monomial_degree(factors),
        total_elementary_factor_offdiagonal_monomials =
            total_elementary_factor_offdiagonal_monomials(factors),
    )
end

function _steinberg_rewrite_count(record, field::Symbol)
    hasproperty(record, field) ||
        throw(ArgumentError("Steinberg rewrite record must include $(field)"))
    value = getproperty(record, field)
    value isa Integer ||
        throw(ArgumentError("Steinberg rewrite record $(field) must be an integer"))
    count = Int(value)
    count >= 0 ||
        throw(ArgumentError("Steinberg rewrite record $(field) must be nonnegative"))
    return count
end

function _steinberg_span_count(span, label::AbstractString)
    span === nothing && return nothing
    hasproperty(span, :start) && hasproperty(span, :stop) ||
        throw(ArgumentError("$(label) span must expose start and stop"))
    start_idx = getproperty(span, :start)
    stop_idx = getproperty(span, :stop)
    start_idx isa Integer && stop_idx isa Integer ||
        throw(ArgumentError("$(label) span bounds must be integers"))
    start_int = Int(start_idx)
    stop_int = Int(stop_idx)
    start_int >= 1 ||
        throw(ArgumentError("$(label) span start must be positive"))
    return stop_int < start_int ? 0 : stop_int - start_int + 1
end

function _normalize_steinberg_rewrite_record(record)
    hasproperty(record, :rule_name) ||
        throw(ArgumentError("Steinberg rewrite record must include rule_name"))
    rule_name = getproperty(record, :rule_name)
    rule_name isa Symbol ||
        throw(ArgumentError("Steinberg rewrite record rule_name must be a Symbol"))
    original_factor_count = _steinberg_rewrite_count(record, :original_factor_count)
    optimized_factor_count = _steinberg_rewrite_count(record, :optimized_factor_count)
    original_span = hasproperty(record, :original_span) ? getproperty(record, :original_span) : nothing
    optimized_span = hasproperty(record, :optimized_span) ? getproperty(record, :optimized_span) : nothing
    metadata = hasproperty(record, :metadata) ? getproperty(record, :metadata) : (;)

    original_span_count = _steinberg_span_count(original_span, "original")
    optimized_span_count = _steinberg_span_count(optimized_span, "optimized")
    original_span_count === nothing || original_span_count == original_factor_count ||
        throw(ArgumentError("Steinberg rewrite original span length must match original_factor_count"))
    optimized_span_count === nothing || optimized_span_count == optimized_factor_count ||
        throw(ArgumentError("Steinberg rewrite optimized span length must match optimized_factor_count"))

    return (;
        rule_name,
        original_factor_count,
        optimized_factor_count,
        original_span,
        optimized_span,
        metadata,
    )
end

function _normalize_steinberg_rewrite_records(records)
    return [_normalize_steinberg_rewrite_record(record) for record in records]
end

function _steinberg_rewrite_log_delta(records)::Int
    return sum(record.optimized_factor_count - record.original_factor_count for record in records; init = 0)
end

function _steinberg_comparison_summary(
    original_factors,
    optimized_factors,
    applied_rewrites,
    original_product,
    optimized_product,
    verification_status::Bool,
)
    original_factor_count = length(original_factors)
    optimized_factor_count = length(optimized_factors)
    return (;
        original_factor_count,
        optimized_factor_count,
        factor_count_delta = optimized_factor_count - original_factor_count,
        original_metrics = _steinberg_metric_summary(original_factors),
        optimized_metrics = _steinberg_metric_summary(optimized_factors),
        applied_rewrites = copy(applied_rewrites),
        original_product,
        optimized_product,
        products_equal = original_product == optimized_product,
        verification_status,
    )
end

function _steinberg_summary_core_status(
    original_factors,
    optimized_factors,
    applied_rewrites,
    original_product,
    optimized_product,
)::Bool
    factor_count_delta = length(optimized_factors) - length(original_factors)
    return original_product == optimized_product &&
           _steinberg_rewrite_log_delta(applied_rewrites) == factor_count_delta
end

function _steinberg_optimization_core_verification(certificate)
    original_context = _steinberg_sequence_context(certificate.original_factors, "original")
    optimized_context =
        _steinberg_sequence_context(certificate.optimized_factors, "optimized", original_context)
    applied_rewrites = _normalize_steinberg_rewrite_records(certificate.applied_rewrites)

    replayed_original_product =
        _steinberg_factor_product(original_context.factors, original_context.ring, original_context.n)
    replayed_optimized_product =
        _steinberg_factor_product(optimized_context.factors, original_context.ring, original_context.n)

    original_product_replay_ok = replayed_original_product == certificate.original_product
    optimized_product_replay_ok = replayed_optimized_product == certificate.optimized_product
    products_equal = certificate.original_product == certificate.optimized_product
    rewrite_log_delta_ok =
        _steinberg_rewrite_log_delta(applied_rewrites) ==
        length(optimized_context.factors) - length(original_context.factors)
    summary_core_status =
        original_product_replay_ok && optimized_product_replay_ok && products_equal && rewrite_log_delta_ok

    expected_summary = _steinberg_comparison_summary(
        original_context.factors,
        optimized_context.factors,
        applied_rewrites,
        certificate.original_product,
        certificate.optimized_product,
        summary_core_status,
    )
    comparison_summary_ok = certificate.comparison_summary == expected_summary
    overall_ok = summary_core_status && comparison_summary_ok

    return (;
        original_sequence_ok = true,
        optimized_sequence_ok = true,
        original_product_replay_ok,
        optimized_product_replay_ok,
        products_equal,
        rewrite_log_delta_ok,
        comparison_summary_ok,
        overall_ok,
    )
end

function _steinberg_optimization_verification(certificate)
    core = _steinberg_optimization_core_verification(certificate)
    stored_verification_ok = certificate.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_ok && stored_verification_ok,
    ))
end

function _steinberg_optimization_certificate(
    original_factors,
    optimized_factors,
    applied_rewrite_metadata = (),
)
    original_context = _steinberg_sequence_context(original_factors, "original")
    optimized_context =
        _steinberg_sequence_context(optimized_factors, "optimized", original_context)
    applied_rewrites = _normalize_steinberg_rewrite_records(applied_rewrite_metadata)
    original_product =
        _steinberg_factor_product(original_context.factors, original_context.ring, original_context.n)
    optimized_product =
        _steinberg_factor_product(optimized_context.factors, original_context.ring, original_context.n)
    summary_status = _steinberg_summary_core_status(
        original_context.factors,
        optimized_context.factors,
        applied_rewrites,
        original_product,
        optimized_product,
    )
    comparison_summary = _steinberg_comparison_summary(
        original_context.factors,
        optimized_context.factors,
        applied_rewrites,
        original_product,
        optimized_product,
        summary_status,
    )
    provisional = SteinbergOptimizationCertificate(
        original_context.factors,
        optimized_context.factors,
        applied_rewrites,
        comparison_summary,
        original_product,
        optimized_product,
        nothing,
    )
    verification = _steinberg_optimization_core_verification(provisional)
    return SteinbergOptimizationCertificate(
        original_context.factors,
        optimized_context.factors,
        applied_rewrites,
        comparison_summary,
        original_product,
        optimized_product,
        verification,
    )
end

function _verify_steinberg_optimization_certificate(certificate)::Bool
    try
        return _steinberg_optimization_verification(certificate).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
