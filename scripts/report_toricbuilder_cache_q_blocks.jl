module ToricBuilderCacheQBlockStatusReport

using Dates
using Serialization
using Oscar
using Printf
using Suslin

include(joinpath(@__DIR__, "..", "test", "fixtures", "toricbuilder_cache_q_blocks.jl"))

const DEFAULT_EXERCISED_CASE_IDS =
    ("case_001", "case_002", "case_003", "case_004", "case_005", "case_006", "case_010")
const DEFAULT_REPORT_PATH =
    joinpath(@__DIR__, "..", "docs", "audits", "2026-06-24-toricbuilder-cache-q-block-status.md")
const SOURCE_FIXTURE = "test/fixtures/toricbuilder_cache_q_blocks.jl"

function _symbol_text(value)
    text = string(value)
    return startswith(text, ":") ? text[2:end] : text
end

function _size_text(size_tuple)
    return string(size_tuple[1], "x", size_tuple[2])
end

function _runtime_text(value)
    value isa Number && return @sprintf("%.3f", value)
    return _symbol_text(value)
end

function _elapsed_seconds(start_ns::UInt64)
    return round((time_ns() - start_ns) / 1.0e9; digits = 3)
end

const STAGE_NAMES = (
    :determinant_classification,
    :normalization,
    :certificate_construction,
    :verification,
)

function _stage_timing(status::Symbol; elapsed_seconds = :not_run, error_details = "none")
    return (; status, elapsed_seconds, error_details)
end

function _not_run_stage_timings()
    return (;
        determinant_classification = _stage_timing(:not_run; error_details = "not_run"),
        normalization = _stage_timing(:not_run; error_details = "not_run"),
        certificate_construction = _stage_timing(:not_run; error_details = "not_run"),
        verification = _stage_timing(:not_run; error_details = "not_run"),
    )
end

function _stage_timings_from_dict(timings::Dict{Symbol, Any})
    return (;
        determinant_classification = get(
            timings,
            :determinant_classification,
            _stage_timing(:not_run; error_details = "not_run"),
        ),
        normalization = get(timings, :normalization, _stage_timing(:not_run; error_details = "not_run")),
        certificate_construction = get(
            timings,
            :certificate_construction,
            _stage_timing(:not_run; error_details = "not_run"),
        ),
        verification = get(timings, :verification, _stage_timing(:not_run; error_details = "not_run")),
    )
end

function _write_worker_progress(
    progress_path,
    current_stage::Symbol,
    stage_started_at::Float64,
    timings;
    peel_progress = nothing,
    determinant_strategy = :eager,
    correction_side = :not_run,
    determinant_source = :not_run,
)
    progress_path === nothing && return nothing
    tmp = string(progress_path, ".tmp")
    open(tmp, "w") do io
        println(io, "current_stage=", current_stage)
        println(io, "stage_started_at=", stage_started_at)
        println(io, "determinant_strategy=", determinant_strategy)
        println(io, "correction_side=", correction_side)
        println(io, "determinant_source=", determinant_source)
        for stage in STAGE_NAMES
            timing = get(timings, stage, nothing)
            timing === nothing && continue
            println(io, "stage.", stage, ".status=", timing.status)
            println(io, "stage.", stage, ".elapsed_seconds=", timing.elapsed_seconds)
            println(io, "stage.", stage, ".error_details=", replace(timing.error_details, '\n' => "\\n"))
        end
        if peel_progress !== nothing
            println(io, "peel.current_dimension=", peel_progress.current_dimension)
            println(io, "peel.completed_steps=", peel_progress.completed_steps)
            peel_progress.last_completed_dimension !== nothing &&
                println(io, "peel.last_completed_dimension=", peel_progress.last_completed_dimension)
            peel_progress.last_completed_elapsed_seconds !== nothing &&
                println(
                    io,
                    "peel.last_completed_elapsed_seconds=",
                    peel_progress.last_completed_elapsed_seconds,
                )
            peel_progress.last_completed_left_factors !== nothing &&
                println(io, "peel.last_completed_left_factors=", peel_progress.last_completed_left_factors)
            peel_progress.last_completed_right_factors !== nothing &&
                println(io, "peel.last_completed_right_factors=", peel_progress.last_completed_right_factors)
            peel_progress.last_column_nnz !== nothing &&
                println(io, "peel.last_column_nnz=", peel_progress.last_column_nnz)
            peel_progress.max_entry_terms !== nothing &&
                println(io, "peel.max_entry_terms=", peel_progress.max_entry_terms)
        end
    end
    mv(tmp, progress_path; force = true)
    return nothing
end

function _maybe_parse_int(raw)
    raw === nothing && return nothing
    return try
        parse(Int, raw)
    catch err
        err isa ArgumentError || rethrow()
        nothing
    end
end

function _maybe_parse_float(raw)
    raw === nothing && return nothing
    return try
        parse(Float64, raw)
    catch err
        err isa ArgumentError || rethrow()
        nothing
    end
end

function _read_peel_progress(data::Dict{String, String})
    current_dimension = _maybe_parse_int(get(data, "peel.current_dimension", nothing))
    completed_steps = _maybe_parse_int(get(data, "peel.completed_steps", nothing))
    last_completed_dimension = _maybe_parse_int(get(data, "peel.last_completed_dimension", nothing))
    last_completed_elapsed_seconds =
        _maybe_parse_float(get(data, "peel.last_completed_elapsed_seconds", nothing))
    last_completed_left_factors = _maybe_parse_int(get(data, "peel.last_completed_left_factors", nothing))
    last_completed_right_factors =
        _maybe_parse_int(get(data, "peel.last_completed_right_factors", nothing))
    last_column_nnz = _maybe_parse_int(get(data, "peel.last_column_nnz", nothing))
    max_entry_terms = _maybe_parse_int(get(data, "peel.max_entry_terms", nothing))
    current_dimension === nothing && return nothing
    completed_steps === nothing && return nothing
    (last_column_nnz === nothing && max_entry_terms === nothing) && return nothing
    return (;
        current_dimension,
        completed_steps,
        last_completed_dimension,
        last_completed_elapsed_seconds,
        last_completed_left_factors,
        last_completed_right_factors,
        last_column_nnz,
        max_entry_terms,
    )
end

function _peel_progress_text(progress)
    progress === nothing && return nothing
    parts = String[
        "peel progress: current d=",
        string(progress.current_dimension),
        ", completed steps=",
        string(progress.completed_steps),
    ]
    if progress.last_completed_dimension !== nothing
        push!(parts, ", last completed d=", string(progress.last_completed_dimension))
        detail_parts = String[]
        progress.last_completed_elapsed_seconds !== nothing &&
            push!(detail_parts, "elapsed $(_runtime_text(progress.last_completed_elapsed_seconds))s")
        progress.last_completed_left_factors !== nothing &&
            push!(detail_parts, "left factors=$(progress.last_completed_left_factors)")
        progress.last_completed_right_factors !== nothing &&
            push!(detail_parts, "right factors=$(progress.last_completed_right_factors)")
        isempty(detail_parts) || push!(parts, " (", join(detail_parts, ", "), ")")
    end
    progress.last_column_nnz !== nothing &&
        push!(parts, ", last-column nnz=", string(progress.last_column_nnz))
    progress.max_entry_terms !== nothing &&
        push!(parts, ", max entry terms=", string(progress.max_entry_terms))
    return join(parts)
end

function _read_worker_progress(progress_path)
    progress_path === nothing &&
        return (;
            current_stage = :determinant_classification,
            timings = Dict{Symbol, Any}(),
            peel_progress = nothing,
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        )
    isfile(progress_path) ||
        return (;
            current_stage = :determinant_classification,
            timings = Dict{Symbol, Any}(),
            peel_progress = nothing,
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        )
    data = Dict{String, String}()
    for line in eachline(progress_path)
        key, value = split(line, "="; limit = 2)
        data[key] = value
    end
    timings = Dict{Symbol, Any}()
    for stage in STAGE_NAMES
        prefix = string("stage.", stage, ".")
        haskey(data, string(prefix, "status")) || continue
        elapsed_raw = data[string(prefix, "elapsed_seconds")]
        elapsed = elapsed_raw == "not_run" ? :not_run : parse(Float64, elapsed_raw)
        timings[stage] = _stage_timing(
            Symbol(data[string(prefix, "status")]);
            elapsed_seconds = elapsed,
            error_details = replace(get(data, string(prefix, "error_details"), "none"), "\\n" => "\n"),
        )
    end
    return (;
        current_stage = Symbol(get(data, "current_stage", "determinant_classification")),
        stage_started_at = parse(Float64, get(data, "stage_started_at", string(time()))),
        timings,
        peel_progress = _read_peel_progress(data),
        determinant_strategy = Symbol(get(data, "determinant_strategy", "eager")),
        correction_side = Symbol(get(data, "correction_side", "not_run")),
        determinant_source = Symbol(get(data, "determinant_source", "not_run")),
    )
end

function _write_peel_worker_progress(
    progress_path,
    stage_started_at::Float64,
    timings::Dict{Symbol, Any},
    progress,
    ;
    determinant_strategy = :eager,
    correction_side = :not_run,
    determinant_source = :not_run,
)
    return _write_worker_progress(
        progress_path,
        :certificate_construction,
        stage_started_at,
        timings;
        peel_progress = progress,
        determinant_strategy,
        correction_side,
        determinant_source,
    )
end

function _record_stage!(timings::Dict{Symbol, Any}, stage::Symbol, f)
    start_ns = time_ns()
    try
        value = f()
        timings[stage] = _stage_timing(:pass; elapsed_seconds = _elapsed_seconds(start_ns))
        return (; status = :pass, value)
    catch err
        err isa InterruptException && rethrow()
        error_details = sprint(showerror, err)
        status = _staged_argument_error(err) ? :certified_algorithm_boundary : :route_error
        timings[stage] = _stage_timing(status; elapsed_seconds = _elapsed_seconds(start_ns), error_details)
        return (; status, error = err, error_details)
    end
end

function _record_worker_stage!(
    f,
    timings::Dict{Symbol, Any},
    stage::Symbol,
    progress_path;
    determinant_strategy = :eager,
    correction_side = :not_run,
    determinant_source = :not_run,
)
    stage_started_at = time()
    _write_worker_progress(
        progress_path,
        stage,
        stage_started_at,
        timings;
        determinant_strategy,
        correction_side,
        determinant_source,
    )
    result = _record_stage!(timings, stage, f)
    _write_worker_progress(
        progress_path,
        stage,
        stage_started_at,
        timings;
        determinant_strategy,
        correction_side,
        determinant_source,
    )
    return result
end

function _stage_timing_text(timing)
    status = _symbol_text(timing.status)
    timing.elapsed_seconds isa Number && return string(status, " (", _runtime_text(timing.elapsed_seconds), "s)")
    return status
end

function _print_stage_timing_details(io, rows)
    println(io, "## Stage Timing Details")
    println(io)
    println(io, "| Case | Determinant classification | Normalization | Certificate construction | Verification |")
    println(io, "| --- | --- | --- | --- | --- |")
    for row in rows
        timings = row.stage_timings
        println(
            io,
            "| $(row.case_id) | $(_stage_timing_text(timings.determinant_classification)) | " *
            "$(_stage_timing_text(timings.normalization)) | " *
            "$(_stage_timing_text(timings.certificate_construction)) | " *
            "$(_stage_timing_text(timings.verification)) |",
        )
    end
    println(io)
end

function _parse_timeout_seconds(raw::AbstractString)
    value = try
        parse(Float64, raw)
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("--timeout-seconds must be a positive number, got $(repr(raw))"))
    end
    isfinite(value) && value > 0 ||
        throw(ArgumentError("--timeout-seconds must be positive, got $(raw)"))
    return value
end

function _parse_choice_symbol(raw::AbstractString, allowed, option_name::AbstractString)
    value = Symbol(raw)
    value in allowed && return value
    allowed_text = join(string.(":", allowed), " or ")
    throw(ArgumentError("$(option_name) must be $(allowed_text), got $(repr(raw))"))
end

function _normalize_report_route_options(determinant_strategy::Symbol, correction_side)
    determinant_strategy in (:eager, :lazy) ||
        throw(ArgumentError("--determinant-strategy must be :eager or :lazy"))
    if determinant_strategy == :eager
        correction_side === nothing ||
            throw(ArgumentError("--correction-side is supported only with --determinant-strategy=lazy"))
        return (; determinant_strategy, correction_side = :not_run)
    end
    resolved_side = correction_side === nothing ? :row : correction_side
    resolved_side in (:row, :column) ||
        throw(ArgumentError("--correction-side must be :row or :column"))
    return (; determinant_strategy, correction_side = resolved_side)
end

function _staged_argument_error(err)
    err isa ArgumentError || return false
    message = sprint(showerror, err)
    return occursin("staged", message) ||
        occursin("unsupported Laurent GL_n determinant", message) ||
        occursin("unsupported exact unimodular column reduction", message) ||
        occursin("normalization boundary", message)
end

function _public_elementary_status(A)
    try
        factors = elementary_factorization(A)
        return verify_factorization(A, factors) ? :elementary_factorization_pass : :factorization_unverified
    catch err
        err isa InterruptException && rethrow()
        return _staged_argument_error(err) ? :staged_boundary : :route_error
    end
end

function _row_route_metadata(;
    determinant_strategy = :eager,
    correction_side = :not_run,
    determinant_source = :not_run,
)
    return (; determinant_strategy, correction_side, determinant_source)
end

function _pending_row(entry)
    return merge((;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = :not_exercised_in_default_report,
        public_elementary_status = :not_run,
        determinant_class = :not_run,
        determinant = "not_run",
        normalization_status = :not_run,
        gl_certificate_status = :not_run,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = :not_run,
        runtime_seconds = :not_run,
        error_details = "not_run",
        evidence = "Fixture recorded; full Suslin route audit not exercised in the default report run.",
        stage_timings = _not_run_stage_timings(),
    ), _row_route_metadata(;
        determinant_strategy = :not_run,
        correction_side = :not_run,
        determinant_source = :not_run,
    ))
end

function _stage_failure_stage(timings::Dict{Symbol, Any}, status::Symbol)
    for stage in STAGE_NAMES
        timing = get(timings, stage, nothing)
        timing !== nothing && timing.status == status && return stage
    end
    return :unknown_stage
end

function _stage_failure_row(
    entry,
    timings,
    result,
    start_ns::UInt64;
    public_status = :not_run,
    profile = nothing,
    determinant_strategy = :eager,
    correction_side = :not_run,
    determinant_source = :not_run,
)
    determinant_class = profile === nothing ? result.status : profile.classification
    determinant = profile === nothing ? string(result.status) : string(profile.determinant)
    return merge((;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = result.status,
        public_elementary_status = public_status,
        determinant_class,
        determinant,
        normalization_status = get(timings, :normalization, _stage_timing(:not_run)).status,
        gl_certificate_status = get(timings, :certificate_construction, _stage_timing(:not_run)).status,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = 0,
        runtime_seconds = _elapsed_seconds(start_ns),
        error_details = result.error_details,
        evidence = "Bounded certificate route stopped at $(result.status) during $(_stage_failure_stage(timings, result.status)); see Route Error Details.",
        stage_timings = _stage_timings_from_dict(timings),
    ), _row_route_metadata(; determinant_strategy, correction_side, determinant_source))
end

function _lazy_gl_certificate_with_progress(
    A;
    correction_side = :row,
    progress_path = nothing,
    stage_started_at = time(),
    timings = Dict{Symbol, Any}(),
    determinant_source_ref = Ref(:not_reached),
)
    progress_callback = progress -> _write_peel_worker_progress(
        progress_path,
        stage_started_at,
        timings,
        progress;
        determinant_strategy = :lazy,
        correction_side,
        determinant_source = determinant_source_ref[],
    )
    deferred_certificate = Suslin._laurent_determinant_deferred_peel_certificate(
        A;
        progress_callback,
    )
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(deferred_certificate)
    determinant_source_ref[] = metadata.determinant_source
    _write_worker_progress(
        progress_path,
        :certificate_construction,
        stage_started_at,
        timings;
        determinant_strategy = :lazy,
        correction_side,
        determinant_source = determinant_source_ref[],
    )
    return Suslin._laurent_gl_lazy_deferred_correction_certificate(
        metadata;
        correction_side,
        progress_callback,
    )
end

function _exercised_row(entry; determinant_strategy = :eager, correction_side = :not_run)
    start_ns = time_ns()
    A = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    public_status = _public_elementary_status(A)
    stage_timings = Dict{Symbol, Any}()

    if determinant_strategy == :lazy
        determinant_source_ref = Ref(:not_reached)
        certificate_result = _record_stage!(
            stage_timings,
            :certificate_construction,
            () -> _lazy_gl_certificate_with_progress(
                A;
                correction_side,
                timings = stage_timings,
                determinant_source_ref,
            ),
        )
        certificate_result.status == :pass ||
            return _stage_failure_row(
                entry,
                stage_timings,
                certificate_result,
                start_ns;
                public_status,
                determinant_strategy,
                correction_side,
                determinant_source = determinant_source_ref[],
            )
        certificate = certificate_result.value

        verification_result = _record_stage!(
            stage_timings,
            :verification,
            () -> verify_laurent_gl_factorization_certificate(certificate),
        )
        verification_result.status == :pass ||
            return _stage_failure_row(
                entry,
                stage_timings,
                verification_result,
                start_ns;
                public_status,
                determinant_strategy,
                correction_side,
                determinant_source = certificate.determinant_source,
            )
        verified = verification_result.value
        return merge((;
            case_id = entry.id,
            matrix_size = entry.dimensions.matrix,
            sparse_entry_count = entry.sparse_entry_count,
            expected_test_level = entry.expected_test_level,
            route_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            public_elementary_status = public_status,
            determinant_class = :deferred_submatrix,
            determinant = string(certificate.overall_determinant),
            normalization_status = :not_run,
            gl_certificate_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            verified,
            factor_count = length(certificate.elementary_factors),
            decomposed_base_matrix_count = length(certificate.elementary_factors),
            runtime_seconds = _elapsed_seconds(start_ns),
            error_details = "none",
            evidence = "Lazy deferred determinant route exercised with $(correction_side) correction; deferred determinant source is $(certificate.determinant_source).",
            stage_timings = _stage_timings_from_dict(stage_timings),
        ), _row_route_metadata(;
            determinant_strategy = :lazy,
            correction_side,
            determinant_source = certificate.determinant_source,
        ))
    end

    profile_result = _record_stage!(stage_timings, :determinant_classification, () -> classify_laurent_determinant(A))
    profile_result.status == :pass ||
        return _stage_failure_row(
            entry,
            stage_timings,
            profile_result,
            start_ns;
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        )
    profile = profile_result.value

    try
        normalization_result =
            _record_stage!(stage_timings, :normalization, () -> normalize_laurent_gl_matrix(A))
        normalization_result.status == :pass ||
            return _stage_failure_row(
                entry,
                stage_timings,
                normalization_result,
                start_ns;
                profile = profile,
                determinant_strategy = :eager,
                correction_side = :not_run,
                determinant_source = :not_run,
            )
        normalization = normalization_result.value

        certificate_result = _record_stage!(
            stage_timings,
            :certificate_construction,
            () -> laurent_gl_factorization_certificate(A),
        )
        certificate_result.status == :pass ||
            return _stage_failure_row(
                entry,
                stage_timings,
                certificate_result,
                start_ns;
                profile = profile,
                determinant_strategy = :eager,
                correction_side = :not_run,
                determinant_source = :not_run,
            )
        certificate = certificate_result.value

        verification_result = _record_stage!(
            stage_timings,
            :verification,
            () -> verify_laurent_gl_factorization_certificate(certificate),
        )
        verification_result.status == :pass ||
            return _stage_failure_row(
                entry,
                stage_timings,
                verification_result,
                start_ns;
                profile = profile,
                determinant_strategy = :eager,
                correction_side = :not_run,
                determinant_source = :not_run,
            )
        verified = verification_result.value
        return merge((;
            case_id = entry.id,
            matrix_size = entry.dimensions.matrix,
            sparse_entry_count = entry.sparse_entry_count,
            expected_test_level = entry.expected_test_level,
            route_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            public_elementary_status = public_status,
            determinant_class = profile.classification,
            determinant = string(profile.determinant),
            normalization_status = :normalization_pass,
            gl_certificate_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            verified,
            factor_count = length(certificate.core_factors),
            decomposed_base_matrix_count = length(certificate.core_factors),
            runtime_seconds = _elapsed_seconds(start_ns),
            error_details = "none",
            evidence = "normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is $(det(normalization.normalized_matrix)).",
            stage_timings = _stage_timings_from_dict(stage_timings),
        ), _row_route_metadata(;
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        ))
    catch err
        err isa InterruptException && rethrow()
        error_details = sprint(showerror, err)
        return merge((;
            case_id = entry.id,
            matrix_size = entry.dimensions.matrix,
            sparse_entry_count = entry.sparse_entry_count,
            expected_test_level = entry.expected_test_level,
            route_status = _staged_argument_error(err) ? :unsupported_staged : :route_error,
            public_elementary_status = public_status,
            determinant_class = profile.classification,
            determinant = string(profile.determinant),
            normalization_status = :normalization_or_certificate_fail,
            gl_certificate_status = :gl_certificate_fail,
            verified = false,
            factor_count = 0,
            decomposed_base_matrix_count = 0,
            runtime_seconds = _elapsed_seconds(start_ns),
            error_details,
            evidence = "Route probe failed; see Route Error Details.",
            stage_timings = _stage_timings_from_dict(stage_timings),
        ), _row_route_metadata(;
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        ))
    end
end

function _worker_exercised_row(
    case_id::AbstractString,
    progress_path;
    determinant_strategy = :eager,
    correction_side = :not_run,
)
    catalog = ToricBuilderCacheQBlocks.catalog()
    _validate_exercised_case_ids!(Set([String(case_id)]), catalog)
    entry = only(filter(entry -> entry.id == case_id, catalog.cases))

    start_ns = time_ns()
    timings = Dict{Symbol, Any}()
    A = ToricBuilderCacheQBlocks.materialize_matrix(entry)

    if determinant_strategy == :lazy
        determinant_source_ref = Ref(:not_reached)
        certificate_stage_started_at = time()
        _write_worker_progress(
            progress_path,
            :certificate_construction,
            certificate_stage_started_at,
            timings;
            determinant_strategy = :lazy,
            correction_side,
            determinant_source = determinant_source_ref[],
        )
        certificate_result = _record_stage!(
            timings,
            :certificate_construction,
            () -> _lazy_gl_certificate_with_progress(
                A;
                correction_side,
                progress_path,
                stage_started_at = certificate_stage_started_at,
                timings,
                determinant_source_ref,
            ),
        )
        determinant_source =
            certificate_result.status == :pass ? certificate_result.value.determinant_source : determinant_source_ref[]
        _write_worker_progress(
            progress_path,
            :certificate_construction,
            certificate_stage_started_at,
            timings;
            determinant_strategy = :lazy,
            correction_side,
            determinant_source,
        )
        certificate_result.status == :pass ||
            return _stage_failure_row(
                entry,
                timings,
                certificate_result,
                start_ns;
                determinant_strategy = :lazy,
                correction_side,
                determinant_source,
            )
        certificate = certificate_result.value

        verification_result = _record_worker_stage!(
            timings,
            :verification,
            progress_path;
            determinant_strategy = :lazy,
            correction_side,
            determinant_source = certificate.determinant_source,
        ) do
            verify_laurent_gl_factorization_certificate(certificate)
        end
        verification_result.status == :pass ||
            return _stage_failure_row(
                entry,
                timings,
                verification_result,
                start_ns;
                determinant_strategy = :lazy,
                correction_side,
                determinant_source = certificate.determinant_source,
            )
        verified = verification_result.value

        return merge((;
            case_id = entry.id,
            matrix_size = entry.dimensions.matrix,
            sparse_entry_count = entry.sparse_entry_count,
            expected_test_level = entry.expected_test_level,
            route_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            public_elementary_status = :not_run,
            determinant_class = :deferred_submatrix,
            determinant = string(certificate.overall_determinant),
            normalization_status = :not_run,
            gl_certificate_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            verified,
            factor_count = length(certificate.elementary_factors),
            decomposed_base_matrix_count = length(certificate.elementary_factors),
            runtime_seconds = _elapsed_seconds(start_ns),
            error_details = "none",
            evidence = "Bounded lazy certificate route exercised with $(correction_side) correction; deferred determinant source is $(certificate.determinant_source).",
            stage_timings = _stage_timings_from_dict(timings),
        ), _row_route_metadata(;
            determinant_strategy = :lazy,
            correction_side,
            determinant_source = certificate.determinant_source,
        ))
    end

    profile_result = _record_worker_stage!(timings, :determinant_classification, progress_path) do
        classify_laurent_determinant(A)
    end
    if profile_result.status != :pass
        return _stage_failure_row(
            entry,
            timings,
            profile_result,
            start_ns;
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        )
    end
    profile = profile_result.value

    normalization_result = _record_worker_stage!(timings, :normalization, progress_path) do
        normalize_laurent_gl_matrix(A)
    end
    normalization_result.status == :pass ||
        return _stage_failure_row(
            entry,
            timings,
            normalization_result,
            start_ns;
            profile,
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        )
    normalization = normalization_result.value

    certificate_stage_started_at = time()
    _write_worker_progress(progress_path, :certificate_construction, certificate_stage_started_at, timings)
    progress_callback = progress -> _write_peel_worker_progress(
        progress_path,
        certificate_stage_started_at,
        timings,
        progress,
        determinant_strategy = :eager,
        correction_side = :not_run,
        determinant_source = :not_run,
    )
    certificate_result = _record_stage!(
        timings,
        :certificate_construction,
        () -> Suslin._laurent_gl_factorization_certificate_from_normalization(
            A,
            normalization;
            progress_callback,
        ),
    )
    _write_worker_progress(progress_path, :certificate_construction, certificate_stage_started_at, timings)
    certificate_result.status == :pass ||
        return _stage_failure_row(
            entry,
            timings,
            certificate_result,
            start_ns;
            profile,
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        )
    certificate = certificate_result.value

    verification_result = _record_worker_stage!(timings, :verification, progress_path) do
        verify_laurent_gl_factorization_certificate(certificate)
    end
    verification_result.status == :pass ||
        return _stage_failure_row(
            entry,
            timings,
            verification_result,
            start_ns;
            profile,
            determinant_strategy = :eager,
            correction_side = :not_run,
            determinant_source = :not_run,
        )
    verified = verification_result.value

    return merge((;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
        public_elementary_status = :not_run,
        determinant_class = profile.classification,
        determinant = string(profile.determinant),
        normalization_status = :normalization_pass,
        gl_certificate_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
        verified,
        factor_count = length(certificate.core_factors),
        decomposed_base_matrix_count = length(certificate.core_factors),
        runtime_seconds = _elapsed_seconds(start_ns),
        error_details = "none",
        evidence = "Bounded certificate route exercised; normalized determinant is $(det(normalization.normalized_matrix)).",
        stage_timings = _stage_timings_from_dict(timings),
    ), _row_route_metadata(;
        determinant_strategy = :eager,
        correction_side = :not_run,
        determinant_source = :not_run,
    ))
end

function _worker_main(
    case_id::AbstractString,
    progress_path,
    result_path;
    determinant_strategy = :eager,
    correction_side = :not_run,
)
    row = _worker_exercised_row(
        case_id,
        progress_path;
        determinant_strategy,
        correction_side,
    )
    if result_path === nothing
        serialize(stdout, row)
    else
        open(result_path, "w") do io
            serialize(io, row)
        end
    end
    return nothing
end

function _catalog_case_ids(catalog)
    return Set(entry.id for entry in catalog.cases)
end

function _validate_exercised_case_ids!(exercised::Set{String}, catalog)
    known = _catalog_case_ids(catalog)
    unknown = sort!(collect(setdiff(exercised, known)))
    if !isempty(unknown)
        known_text = join(sort!(collect(known)), ", ")
        throw(ArgumentError(
            "unknown --exercise case ID(s): $(join(unknown, ", ")); known case IDs: $(known_text)",
        ))
    end
    return exercised
end

function _timed_out_row(entry, timeout_seconds, runtime_seconds, progress; cleanup_details = nothing)
    timings = copy(progress.timings)
    current_stage = progress.current_stage
    peel_progress_text =
        hasproperty(progress, :peel_progress) ? _peel_progress_text(progress.peel_progress) : nothing
    stage_elapsed =
        hasproperty(progress, :stage_started_at) ? max(round(time() - progress.stage_started_at; digits = 3), 0.0) :
        runtime_seconds
    stage_error_details = "timed out after $(_runtime_text(timeout_seconds)) seconds"
    if cleanup_details !== nothing
        stage_error_details = string(stage_error_details, "; ", cleanup_details)
    end
    if peel_progress_text !== nothing
        stage_error_details = string(stage_error_details, "; ", peel_progress_text)
    end
    row_error_details = string(stage_error_details, " while running ", current_stage)
    evidence = "Bounded worker exceeded the configured process-level budget."
    if cleanup_details !== nothing
        evidence = string(evidence, " Cleanup warning: ", cleanup_details)
    end
    if peel_progress_text !== nothing
        evidence = string(evidence, " ", peel_progress_text)
    end
    timings[current_stage] = _stage_timing(
        :timed_out;
        elapsed_seconds = stage_elapsed,
        error_details = stage_error_details,
    )
    return merge((;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = :timed_out,
        public_elementary_status = :not_run,
        determinant_class = :timed_out,
        determinant = "timed_out",
        normalization_status = current_stage == :normalization ? :timed_out : get(
            timings,
            :normalization,
            _stage_timing(:not_run),
        ).status,
        gl_certificate_status = current_stage in (:certificate_construction, :verification) ? :timed_out : :not_run,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = 0,
        runtime_seconds,
        error_details = row_error_details,
        evidence,
        stage_timings = _stage_timings_from_dict(timings),
    ), _row_route_metadata(;
        determinant_strategy = hasproperty(progress, :determinant_strategy) ? progress.determinant_strategy : :eager,
        correction_side = hasproperty(progress, :correction_side) ? progress.correction_side : :not_run,
        determinant_source = hasproperty(progress, :determinant_source) ? progress.determinant_source : :not_run,
    ))
end

function _worker_command(
    entry,
    progress_path,
    result_path;
    determinant_strategy = :eager,
    correction_side = :not_run,
)
    project_path = dirname(Base.active_project())
    args = String[
        "--bounded-worker=$(entry.id)",
        "--worker-progress=$(progress_path)",
        "--worker-result=$(result_path)",
        "--determinant-strategy=$(_symbol_text(determinant_strategy))",
    ]
    determinant_strategy == :lazy &&
        push!(args, "--correction-side=$(_symbol_text(correction_side))")
    return Cmd(vcat(
        collect(Base.julia_cmd()),
        ["--project=$(project_path)", abspath(@__FILE__)],
        args,
    ))
end

function _worker_route_error_row(
    entry,
    runtime_seconds,
    stderr_text;
    determinant_strategy = :eager,
    correction_side = :not_run,
    progress = nothing,
)
    error_details = isempty(stderr_text) ? "bounded worker exited nonzero" : stderr_text
    if progress === nothing
        current_stage = :determinant_classification
        timings = Dict{Symbol, Any}()
        route_determinant_strategy = determinant_strategy
        route_correction_side = correction_side
        route_determinant_source = determinant_strategy == :lazy ? :not_reached : :not_run
    else
        current_stage = progress.current_stage
        timings = copy(progress.timings)
        route_determinant_strategy = hasproperty(progress, :determinant_strategy) ?
            progress.determinant_strategy : determinant_strategy
        route_correction_side = hasproperty(progress, :correction_side) ?
            progress.correction_side : correction_side
        route_determinant_source = hasproperty(progress, :determinant_source) ?
            progress.determinant_source :
            (route_determinant_strategy == :lazy ? :not_reached : :not_run)
    end
    timings[current_stage] = _stage_timing(:route_error; elapsed_seconds = runtime_seconds, error_details)
    return merge((;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = :route_error,
        public_elementary_status = :not_run,
        determinant_class = :route_error,
        determinant = "route_error",
        normalization_status = :not_run,
        gl_certificate_status = :not_run,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = 0,
        runtime_seconds,
        error_details,
        evidence = "Bounded worker exited before producing a report row.",
        stage_timings = _stage_timings_from_dict(timings),
    ), _row_route_metadata(;
        determinant_strategy = route_determinant_strategy,
        correction_side = route_correction_side,
        determinant_source = route_determinant_source,
    ))
end

function _wait_for_exit_after_kill(proc; grace_seconds = 1.0, poll_seconds = 0.05)
    deadline = time() + grace_seconds
    while process_running(proc)
        remaining = deadline - time()
        remaining <= 0 && break
        sleep(min(poll_seconds, remaining))
    end
    return !process_running(proc)
end

function _deserialize_worker_row(stdout_path)
    row = open(deserialize, stdout_path)
    hasproperty(row, :route_status) && hasproperty(row, :stage_timings) && return row
    throw(ArgumentError("bounded worker produced invalid stdout payload of type $(typeof(row))"))
end

function _bounded_exercised_row(
    entry,
    timeout_seconds::Float64;
    determinant_strategy = :eager,
    correction_side = :not_run,
)
    stdout_path = tempname()
    stderr_path = tempname()
    progress_path = tempname()
    result_path = tempname()
    start_time = time()
    proc = run(
        pipeline(
            _worker_command(
                entry,
                progress_path,
                result_path;
                determinant_strategy,
                correction_side,
            );
            stdout = stdout_path,
            stderr = stderr_path,
        ),
        wait = false,
    )

    try
        while process_running(proc)
            runtime_seconds = round(time() - start_time; digits = 3)
            if runtime_seconds >= timeout_seconds
                try
                    kill(proc)
                catch
                end
                exited_after_kill = _wait_for_exit_after_kill(proc)
                progress = _read_worker_progress(progress_path)
                cleanup_details = exited_after_kill ? nothing : "worker did not exit after kill grace period"
                return _timed_out_row(
                    entry,
                    timeout_seconds,
                    runtime_seconds,
                    progress;
                    cleanup_details,
                )
            end
            sleep(0.05)
        end

        wait(proc)
        runtime_seconds = round(time() - start_time; digits = 3)
        if success(proc)
            try
                return _deserialize_worker_row(result_path)
            catch err
                err isa InterruptException && rethrow()
                stderr_text = isfile(stderr_path) ? read(stderr_path, String) : ""
                deserialize_error = "deserialize failed: $(sprint(showerror, err))"
                combined_details =
                    isempty(stderr_text) ? deserialize_error : string(stderr_text, "\n", deserialize_error)
                return _worker_route_error_row(
                    entry,
                    runtime_seconds,
                    combined_details;
                    determinant_strategy,
                    correction_side,
                    progress = _read_worker_progress(progress_path),
                )
            end
        end
        stderr_text = isfile(stderr_path) ? read(stderr_path, String) : ""
        return _worker_route_error_row(
            entry,
            runtime_seconds,
            stderr_text;
            determinant_strategy,
            correction_side,
            progress = _read_worker_progress(progress_path),
        )
    finally
        for path in (stdout_path, stderr_path, progress_path, string(progress_path, ".tmp"), result_path)
            isfile(path) && rm(path; force = true)
        end
    end
end

function build_report(;
    exercised_case_ids = DEFAULT_EXERCISED_CASE_IDS,
    generated_on = Dates.today(),
    timeout_seconds = nothing,
    determinant_strategy = :eager,
    correction_side = nothing,
)
    options = _normalize_report_route_options(determinant_strategy, correction_side)
    exercised = Set(string.(exercised_case_ids))
    catalog = ToricBuilderCacheQBlocks.catalog()
    _validate_exercised_case_ids!(exercised, catalog)
    rows = [
        entry.id in exercised ?
        (
            timeout_seconds === nothing ?
            _exercised_row(
                entry;
                determinant_strategy = options.determinant_strategy,
                correction_side = options.correction_side,
            ) :
            _bounded_exercised_row(
                entry,
                timeout_seconds;
                determinant_strategy = options.determinant_strategy,
                correction_side = options.correction_side,
            )
        ) :
        _pending_row(entry)
        for entry in catalog.cases
    ]
    return (;
        title = "ToricBuilder Cache Q-Block Status Report",
        generated_on,
        source_fixture = SOURCE_FIXTURE,
        exercised_case_ids = Tuple(sort!(collect(exercised))),
        timeout_seconds,
        rows,
    )
end

function render_markdown(report)
    io = IOBuffer()
    println(io, "# $(report.title)")
    println(io)
    println(io, "Date: $(report.generated_on)")
    println(io, "Source fixture: `$(report.source_fixture)`")
    println(io)
    println(io, "This report records the current Suslin status for the checked-in ToricBuilder cache Q-block fixtures.")
    println(io, "The default run fully exercises `$(join(report.exercised_case_ids, "`, `"))`; every other case is saved fixture data marked `not_exercised_in_default_report`.")
    println(io)
    println(io, "## Summary")
    println(io)
    exercised_count = count(row -> row.public_elementary_status != :not_run, report.rows)
    pass_count = count(row -> row.route_status == :gl_certificate_pass, report.rows)
    error_count = count(row -> row.route_status == :route_error, report.rows)
    not_exercised_count = count(row -> row.route_status == :not_exercised_in_default_report, report.rows)
    println(io, "- Total checked-in Q-block cases: $(length(report.rows))")
    println(io, "- Fully exercised cases: $(exercised_count)")
    println(io, "- GL certificate passes: $(pass_count)")
    println(io, "- Route errors: $(error_count)")
    println(io, "- Not exercised in default report: $(not_exercised_count)")
    println(io)
    println(io, "## Case Table")
    println(io)
    println(io, "| Case | Matrix size | Sparse nnz | Test level | Route status | Public elementary status | Determinant class | Decomposed base matrices | Runtime seconds |")
    println(io, "| --- | ---: | ---: | --- | --- | --- | --- | ---: | ---: |")
    for row in report.rows
        println(
            io,
            "| $(row.case_id) | $(_size_text(row.matrix_size)) | $(row.sparse_entry_count) | " *
            "$(_symbol_text(row.expected_test_level)) | $(_symbol_text(row.route_status)) | " *
            "$(_symbol_text(row.public_elementary_status)) | $(_symbol_text(row.determinant_class)) | " *
            "$(_symbol_text(row.decomposed_base_matrix_count)) | $(_runtime_text(row.runtime_seconds)) |",
        )
    end
    println(io)
    _print_determinant_route_metadata(io, report.rows)
    println(io, "## Exercised Evidence")
    println(io)
    for row in report.rows
        row.public_elementary_status != :not_run || continue
        println(
            io,
            "- `$(row.case_id)`: determinant `$(row.determinant)`, route `$(row.route_status)`, public `$(row.public_elementary_status)`, decomposed base matrices `$(row.decomposed_base_matrix_count)`, verified `$(row.verified)`, determinant_strategy `$(_symbol_text(row.determinant_strategy))`, correction_side `$(_symbol_text(row.correction_side))`, determinant_source `$(_symbol_text(row.determinant_source))`. $(row.evidence)",
        )
    end
    println(io)
    _print_stage_timing_details(io, report.rows)
    error_rows = filter(row -> row.error_details != "none" && row.error_details != "not_run", report.rows)
    if !isempty(error_rows)
        println(io, "## Route Error Details")
        println(io)
        for row in error_rows
            println(io, "- `$(row.case_id)`: $(row.error_details)")
        end
        println(io)
    end
    println(io, "## Not Exercised Boundary")
    println(io)
    println(io, "The non-exercised rows are real stored matrices, not placeholders. They are intentionally not routed through the full factorization stack in the default report because large Laurent GL certification can be much slower than fixture/schema validation.")
    println(io, "Use `--exercise=case_001,case_003,...` to expand the set of cases that the report probes.")
    println(io)
    println(io, "## Reproduction")
    println(io)
    println(io, "```text")
    println(io, "julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl")
    println(io, "julia --project=. test/internal/toricbuilder_cache_status_report.jl")
    println(io, "```")
    return String(take!(io))
end

function _print_determinant_route_metadata(io, rows)
    println(io, "## Determinant Route Metadata")
    println(io)
    println(io, "| Case | determinant_strategy | correction_side | determinant_source |")
    println(io, "| --- | --- | --- | --- |")
    for row in rows
        println(
            io,
            "| $(row.case_id) | $(_symbol_text(row.determinant_strategy)) | $(_symbol_text(row.correction_side)) | $(_symbol_text(row.determinant_source)) |",
        )
    end
    println(io)
end

function write_report(path::AbstractString = DEFAULT_REPORT_PATH; report = build_report())
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, render_markdown(report))
    end
    return path
end

function _parse_args(args)
    output = DEFAULT_REPORT_PATH
    exercised = collect(DEFAULT_EXERCISED_CASE_IDS)
    timeout_seconds = nothing
    determinant_strategy = :eager
    correction_side = nothing

    for arg in args
        if startswith(arg, "--output=")
            output = arg[length("--output=")+1:end]
        elseif startswith(arg, "--exercise=")
            raw = split(arg[length("--exercise=")+1:end], ",")
            exercised = [strip(case_id) for case_id in raw if !isempty(strip(case_id))]
        elseif startswith(arg, "--timeout-seconds=")
            timeout_seconds = _parse_timeout_seconds(arg[length("--timeout-seconds=")+1:end])
        elseif startswith(arg, "--determinant-strategy=")
            determinant_strategy = _parse_choice_symbol(
                arg[length("--determinant-strategy=")+1:end],
                (:eager, :lazy),
                "--determinant-strategy",
            )
        elseif startswith(arg, "--correction-side=")
            correction_side = _parse_choice_symbol(
                arg[length("--correction-side=")+1:end],
                (:row, :column),
                "--correction-side",
            )
        else
            throw(ArgumentError("unsupported argument: $(arg)"))
        end
    end

    options = _normalize_report_route_options(determinant_strategy, correction_side)
    return (;
        output,
        exercised,
        timeout_seconds,
        determinant_strategy = options.determinant_strategy,
        correction_side = options.correction_side,
    )
end

function _worker_arg_value(args, prefix::AbstractString)
    matches = [arg[length(prefix)+1:end] for arg in args if startswith(arg, prefix)]
    return isempty(matches) ? nothing : only(matches)
end

function main(args = ARGS)
    worker_case = _worker_arg_value(args, "--bounded-worker=")
    if worker_case !== nothing
        progress_path = _worker_arg_value(args, "--worker-progress=")
        result_path = _worker_arg_value(args, "--worker-result=")
        determinant_strategy = _parse_choice_symbol(
            something(_worker_arg_value(args, "--determinant-strategy="), "eager"),
            (:eager, :lazy),
            "--determinant-strategy",
        )
        correction_side_raw = _worker_arg_value(args, "--correction-side=")
        correction_side =
            correction_side_raw === nothing ? nothing : _parse_choice_symbol(
                correction_side_raw,
                (:row, :column),
                "--correction-side",
            )
        options = _normalize_report_route_options(determinant_strategy, correction_side)
        _worker_main(
            worker_case,
            progress_path,
            result_path;
            determinant_strategy = options.determinant_strategy,
            correction_side = options.correction_side,
        )
        return nothing
    end

    options = _parse_args(args)
    report = build_report(;
        exercised_case_ids = options.exercised,
        timeout_seconds = options.timeout_seconds,
        determinant_strategy = options.determinant_strategy,
        correction_side = options.correction_side == :not_run ? nothing : options.correction_side,
    )
    path = write_report(options.output; report)
    println("wrote ToricBuilder cache Q-block status report to $(path)")
    return path
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    ToricBuilderCacheQBlockStatusReport.main()
end
