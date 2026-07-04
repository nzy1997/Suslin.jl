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
    (R isa MPolyRing || R isa PolyRing) ||
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

function _steinberg_rewrite_log_counts_within_sequence_lengths(
    records,
    original_factor_count::Int,
    optimized_factor_count::Int,
)::Bool
    total_original_factor_count =
        sum(record.original_factor_count for record in records; init = 0)
    total_optimized_factor_count =
        sum(record.optimized_factor_count for record in records; init = 0)
    return total_original_factor_count <= original_factor_count &&
           total_optimized_factor_count <= optimized_factor_count
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
    rewrite_log_counts_ok = _steinberg_rewrite_log_counts_within_sequence_lengths(
        applied_rewrites,
        length(original_factors),
        length(optimized_factors),
    )
    return original_product == optimized_product &&
           _steinberg_rewrite_log_delta(applied_rewrites) == factor_count_delta &&
           rewrite_log_counts_ok
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
    rewrite_log_counts_ok = _steinberg_rewrite_log_counts_within_sequence_lengths(
        applied_rewrites,
        length(original_context.factors),
        length(optimized_context.factors),
    )
    summary_core_status =
        original_product_replay_ok &&
        optimized_product_replay_ok &&
        products_equal &&
        rewrite_log_delta_ok &&
        rewrite_log_counts_ok

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
        rewrite_log_counts_ok,
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

function _steinberg_same_elementary_position(left, right)::Bool
    return left.kind == :elementary &&
           right.kind == :elementary &&
           left.row == right.row &&
           left.col == right.col
end

function _steinberg_adjacent_rewrite_record(
    rule_name::Symbol,
    original_start::Int,
    original_stop::Int,
    optimized_start::Int,
    optimized_stop::Int,
    metadata,
)
    return (;
        rule_name,
        original_factor_count = original_stop - original_start + 1,
        optimized_factor_count = optimized_stop < optimized_start ? 0 :
            optimized_stop - optimized_start + 1,
        original_span = (; start = original_start, stop = original_stop),
        optimized_span = (; start = optimized_start, stop = optimized_stop),
        metadata,
    )
end

function _steinberg_commutator_inverse_tail_matches(first, second, third, fourth)::Bool
    return third.row == first.row &&
           third.col == first.col &&
           third.coefficient == -first.coefficient &&
           fourth.row == second.row &&
           fourth.col == second.col &&
           fourth.coefficient == -second.coefficient
end

function _steinberg_commutator_forward_candidate(first, second)
    i = first.row
    j = first.col
    l = second.col
    first.col == second.row || return nothing
    i != l || return nothing

    return (;
        rule_name = :commutator_forward,
        replacement_records = (;
            kind = :elementary,
            n = first.n,
            ring = first.ring,
            row = i,
            col = l,
            coefficient = first.coefficient * second.coefficient,
        ),
        metadata = (;
            indices = (; i, j, l),
            a = first.coefficient,
            b = second.coefficient,
        ),
    )
end

function _steinberg_commutator_reverse_candidate(first, second)
    l = second.row
    i = first.row
    j = first.col
    second.col == i || return nothing
    j != l || return nothing

    return (;
        rule_name = :commutator_reverse,
        replacement_records = (;
            kind = :elementary,
            n = first.n,
            ring = first.ring,
            row = l,
            col = j,
            coefficient = -(first.coefficient * second.coefficient),
        ),
        metadata = (;
            indices = (; l, i, j),
            a = first.coefficient,
            b = second.coefficient,
        ),
    )
end

function _steinberg_commutator_disjoint_candidate(first, second)
    i = first.row
    j = first.col
    l = second.row
    p = second.col
    i != p || return nothing
    j != l || return nothing

    return (;
        rule_name = :disjoint_commutator_identity,
        replacement_records = (),
        metadata = (;
            indices = (; i, j, l, p),
            a = first.coefficient,
            b = second.coefficient,
        ),
    )
end

function _steinberg_commutator_window_candidate(window_records)
    all(record -> record.kind == :elementary, window_records) || return nothing
    first, second, third, fourth = window_records
    _steinberg_commutator_inverse_tail_matches(first, second, third, fourth) ||
        return nothing

    candidate = _steinberg_commutator_forward_candidate(first, second)
    candidate === nothing || return candidate
    candidate = _steinberg_commutator_reverse_candidate(first, second)
    candidate === nothing || return candidate
    return _steinberg_commutator_disjoint_candidate(first, second)
end

function _steinberg_commutator_replacement_factors(candidate)
    replacement_records = candidate.replacement_records
    records = hasproperty(replacement_records, :kind) ?
        (replacement_records,) :
        replacement_records
    return [_elementary_factor_record_matrix(record) for record in records]
end

function _steinberg_commutator_local_products_equal(
    original_window_factors,
    replacement_factors,
    R,
    n::Int,
)::Bool
    return _steinberg_factor_product(original_window_factors, R, n) ==
           _steinberg_factor_product(replacement_factors, R, n)
end

function _steinberg_commutator_rewrite_optimization_certificate(factors)
    original_context = _steinberg_sequence_context(factors, "original")
    optimized_factors = Any[]
    applied_rewrites = Any[]
    records = original_context.records

    index = 1
    while index <= length(records)
        if index + 3 <= length(records)
            window_records = records[index:(index + 3)]
            candidate = _steinberg_commutator_window_candidate(window_records)
            if candidate !== nothing
                original_window_factors = original_context.factors[index:(index + 3)]
                replacement_factors = _steinberg_commutator_replacement_factors(candidate)
                local_products_equal = _steinberg_commutator_local_products_equal(
                    original_window_factors,
                    replacement_factors,
                    original_context.ring,
                    original_context.n,
                )

                if local_products_equal
                    optimized_start = length(optimized_factors) + 1
                    append!(optimized_factors, replacement_factors)
                    optimized_stop = optimized_start + length(replacement_factors) - 1
                    push!(
                        applied_rewrites,
                        _steinberg_adjacent_rewrite_record(
                            candidate.rule_name,
                            index,
                            index + 3,
                            optimized_start,
                            optimized_stop,
                            merge(candidate.metadata, (; local_products_equal,)),
                        ),
                    )
                    index += 4
                    continue
                end
            end
        end

        push!(optimized_factors, original_context.factors[index])
        index += 1
    end

    return _steinberg_optimization_certificate(
        original_context.factors,
        optimized_factors,
        applied_rewrites,
    )
end

function _steinberg_adjacent_rewrite_optimization_certificate(factors)
    original_context = _steinberg_sequence_context(factors, "original")
    optimized_factors = Any[]
    applied_rewrites = Any[]
    records = original_context.records

    index = 1
    while index <= length(records)
        record = records[index]

        if record.kind == :identity
            optimized_start = length(optimized_factors) + 1
            push!(
                applied_rewrites,
                _steinberg_adjacent_rewrite_record(
                    :identity_removal,
                    index,
                    index,
                    optimized_start,
                    optimized_start - 1,
                    (; kind = :identity),
                ),
            )
            index += 1
            continue
        end

        run_stop = index
        while run_stop < length(records) &&
              _steinberg_same_elementary_position(record, records[run_stop + 1])
            run_stop += 1
        end

        if run_stop > index
            merged_coefficient = zero(record.ring)
            for run_index in index:run_stop
                merged_coefficient += records[run_index].coefficient
            end

            optimized_start = length(optimized_factors) + 1
            optimized_stop = optimized_start - 1
            rule_name = :same_position_merge
            if !iszero(merged_coefficient)
                merged_record = (;
                    kind = :elementary,
                    n = record.n,
                    ring = record.ring,
                    row = record.row,
                    col = record.col,
                    coefficient = merged_coefficient,
                )
                push!(optimized_factors, _elementary_factor_record_matrix(merged_record))
                optimized_stop = optimized_start
            elseif run_stop == index + 1
                rule_name = :inverse_cancellation
            end

            push!(
                applied_rewrites,
                _steinberg_adjacent_rewrite_record(
                    rule_name,
                    index,
                    run_stop,
                    optimized_start,
                    optimized_stop,
                    (;
                        row = record.row,
                        col = record.col,
                        original_indices = (; start = index, stop = run_stop),
                        merged_coefficient,
                    ),
                ),
            )
            index = run_stop + 1
            continue
        end

        push!(optimized_factors, _elementary_factor_record_matrix(record))
        index += 1
    end

    return _steinberg_optimization_certificate(
        original_context.factors,
        optimized_factors,
        applied_rewrites,
    )
end
