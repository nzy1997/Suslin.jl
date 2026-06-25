module ToricBuilderCacheQBlockStatusReport

using Dates
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

function _staged_argument_error(err)
    err isa ArgumentError || return false
    message = sprint(showerror, err)
    return occursin("staged", message) ||
        occursin("unsupported Laurent GL_n determinant", message) ||
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

function _pending_row(entry)
    return (;
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
    )
end

function _exercised_row(entry)
    start_ns = time_ns()
    A = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    profile = classify_laurent_determinant(A)
    public_status = _public_elementary_status(A)

    try
        normalization = normalize_laurent_gl_matrix(A)
        certificate = laurent_gl_factorization_certificate(A)
        verified = verify_laurent_gl_factorization_certificate(certificate)
        return (;
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
        )
    catch err
        err isa InterruptException && rethrow()
        error_details = sprint(showerror, err)
        return (;
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
        )
    end
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

function build_report(;
    exercised_case_ids = DEFAULT_EXERCISED_CASE_IDS,
    generated_on = Dates.today(),
    timeout_seconds = nothing,
)
    exercised = Set(string.(exercised_case_ids))
    catalog = ToricBuilderCacheQBlocks.catalog()
    _validate_exercised_case_ids!(exercised, catalog)
    rows = [
        entry.id in exercised ? _exercised_row(entry) : _pending_row(entry)
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
    println(io, "## Exercised Evidence")
    println(io)
    for row in report.rows
        row.public_elementary_status != :not_run || continue
        println(io, "- `$(row.case_id)`: determinant `$(row.determinant)`, route `$(row.route_status)`, public `$(row.public_elementary_status)`, decomposed base matrices `$(row.decomposed_base_matrix_count)`, verified `$(row.verified)`. $(row.evidence)")
    end
    println(io)
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

    for arg in args
        if startswith(arg, "--output=")
            output = arg[length("--output=")+1:end]
        elseif startswith(arg, "--exercise=")
            raw = split(arg[length("--exercise=")+1:end], ",")
            exercised = [strip(case_id) for case_id in raw if !isempty(strip(case_id))]
        elseif startswith(arg, "--timeout-seconds=")
            timeout_seconds = _parse_timeout_seconds(arg[length("--timeout-seconds=")+1:end])
        else
            throw(ArgumentError("unsupported argument: $(arg)"))
        end
    end

    return (; output, exercised, timeout_seconds)
end

function main(args = ARGS)
    options = _parse_args(args)
    report = build_report(;
        exercised_case_ids = options.exercised,
        timeout_seconds = options.timeout_seconds,
    )
    path = write_report(options.output; report)
    println("wrote ToricBuilder cache Q-block status report to $(path)")
    return path
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    ToricBuilderCacheQBlockStatusReport.main()
end
