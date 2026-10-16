using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "test", "fixtures", "toricbuilder_case008_d14_column_boundary.jl"))

const CASE008_D14_LTP_SELECTED_ENTRY = 9
const CASE008_D14_LTP_SELECTED_GENERATOR_NAME = "v"
const CASE008_D14_LTP_EXPECTED_CONVERSION_FACTORS = 17
const CASE008_D14_LTP_EXPECTED_SOURCE_TERMS = 50311
const CASE008_D14_LTP_EXPECTED_POLYNOMIAL_TERMS = 2774551
const CASE008_D14_LTP_EXPECTED_POLYNOMIALIZATION_POWER = 93
const CASE008_D14_LTP_REPORT_STATUS = :laurent_to_polynomial_conversion_report
const CASE008_D14_LTP_REPORT_FIELDS = (
    :case_id,
    :source_dimension,
    :ring_generators,
    :source_boundary,
    :selected_entry,
    :selected_generator,
    :conversion_status,
    :conversion_replay,
    :conversion_factors,
    :source_terms,
    :polynomial_terms,
    :polynomialization_power,
    :polynomial_unimodular,
    :polynomial_witness_replay,
    :polynomial_witness_terms,
    :polynomial_witness_provenance,
    :negative_controls,
    :negative_controls_pass,
    :elapsed_seconds,
    :timeout_seconds,
    :final_status,
    :status,
)
const CASE008_D14_LTP_REPORT_CACHE = IdDict{Any, Any}()

function _case008_d14_ltp_cached_report()
    isempty(CASE008_D14_LTP_REPORT_CACHE) && return nothing
    return first(values(CASE008_D14_LTP_REPORT_CACHE))
end

function _case008_d14_ltp_default_fixture()
    cached = _case008_d14_ltp_cached_report()
    cached !== nothing && return cached.fixture
    return ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
end

function _case008_d14_ltp_ring_generator_names(R)
    return Tuple(string.(gens(R)))
end

function _case008_d14_ltp_selected_generator(R)
    generators = tuple(gens(R)...)
    selected = findfirst(generator -> string(generator) == CASE008_D14_LTP_SELECTED_GENERATOR_NAME, generators)
    selected === nothing &&
        throw(ArgumentError("case_008 d14 fixture does not expose selected generator $(CASE008_D14_LTP_SELECTED_GENERATOR_NAME)"))
    return generators[selected]
end

function _case008_d14_ltp_tamper(certificate; kwargs...)
    fields = propertynames(certificate)
    values = map(field -> get(kwargs, field, getproperty(certificate, field)), fields)
    return Suslin.LaurentToPolynomialColumnCertificate(values...)
end

function _case008_d14_ltp_replace_tuple_entry(values::Tuple, idx::Int, value)
    replacement = collect(values)
    replacement[idx] = value
    return tuple(replacement...)
end

function _case008_d14_ltp_modified_preprocessing_status(certificate)
    R = certificate.ring
    n = length(certificate.original_column)
    modified = collect(certificate.conversion_factors)
    idx = findfirst(
        factor -> Suslin._canonical_elementary_factor_record(factor).kind == :elementary,
        modified,
    )
    idx === nothing && return :missing_preprocessing_factor
    record = Suslin._canonical_elementary_factor_record(modified[idx])
    modified[idx] = elementary_matrix(
        n,
        record.row,
        record.col,
        record.coefficient + first(gens(R)),
        R,
    )
    product = Suslin._laurent_to_polynomial_factor_sequence_product(modified, R, n)
    return Suslin._validate_laurent_to_polynomial_certificate(
        _case008_d14_ltp_tamper(
            certificate;
            conversion_factors = modified,
            conversion_product = product,
        ),
    )
end

function _case008_d14_ltp_stale_intermediate_status(certificate)
    R = certificate.ring
    stale = _case008_d14_ltp_replace_tuple_entry(
        certificate.intermediate_laurent_column,
        1,
        certificate.intermediate_laurent_column[1] + one(R),
    )
    return Suslin._validate_laurent_to_polynomial_certificate(
        _case008_d14_ltp_tamper(certificate; intermediate_laurent_column = stale),
    )
end

function _case008_d14_ltp_noninvertible_inverse_status(certificate)
    R = certificate.ring
    bad_inverse = _case008_d14_ltp_replace_tuple_entry(
        certificate.inverse_lift,
        1,
        merge(certificate.inverse_lift[1], (; laurent_value = first(gens(R)) + one(R))),
    )
    return Suslin._validate_laurent_to_polynomial_certificate(
        _case008_d14_ltp_tamper(certificate; inverse_lift = bad_inverse),
    )
end

function _case008_d14_ltp_copied_polynomial_column_status(certificate)
    P = certificate.polynomial_ring
    copied = _case008_d14_ltp_replace_tuple_entry(
        certificate.polynomial_column,
        1,
        certificate.polynomial_column[1] + one(P),
    )
    return Suslin._validate_laurent_to_polynomial_certificate(
        _case008_d14_ltp_tamper(certificate; polynomial_column = copied),
    )
end

function _case008_d14_ltp_modified_witness_status(certificate)
    P = certificate.polynomial_ring
    idx = findfirst(!iszero, certificate.polynomial_column)
    idx === nothing && return :missing_nonzero_polynomial_entry
    witness = collect(certificate.polynomial_witness)
    witness[idx] += one(P)
    return Suslin._validate_laurent_to_polynomial_certificate(
        _case008_d14_ltp_tamper(certificate; polynomial_witness = tuple(witness...)),
    )
end

function _case008_d14_ltp_copied_witness_status(certificate)
    P = certificate.polynomial_ring
    witness = [zero(P) for _ in certificate.polynomial_witness]
    witness[1] = one(P)
    return Suslin._validate_laurent_to_polynomial_certificate(
        _case008_d14_ltp_tamper(certificate; polynomial_witness = tuple(witness...)),
    )
end

function _case008_d14_ltp_nonunimodular_source_status(fixture)
    corrupted = ToricBuilderCase008D14ColumnBoundary.non_unimodular_negative_control(fixture)
    selected_generator = _case008_d14_ltp_selected_generator(corrupted.ring)
    try
        noether = Suslin._laurent_noether_certificate(
            corrupted.failing_column,
            CASE008_D14_LTP_SELECTED_ENTRY,
            selected_generator,
        )
        Suslin._laurent_to_polynomial_certificate(
            corrupted.failing_column,
            noether,
            CASE008_D14_LTP_SELECTED_ENTRY,
            selected_generator,
        )
        return :accepted
    catch err
        err isa InterruptException && rethrow()
        return :rejected
    end
end

function _case008_d14_ltp_negative_controls(fixture, certificate)
    controls = (;
        nonunimodular_source = _case008_d14_ltp_nonunimodular_source_status(fixture),
        modified_preprocessing_factor =
            _case008_d14_ltp_modified_preprocessing_status(certificate),
        stale_intermediate_column =
            _case008_d14_ltp_stale_intermediate_status(certificate),
        modified_polynomialization_power =
            Suslin._validate_laurent_to_polynomial_certificate(
                _case008_d14_ltp_tamper(
                    certificate;
                    polynomialization_power = certificate.polynomialization_power + 1,
                ),
            ),
        noninvertible_inverse_map =
            _case008_d14_ltp_noninvertible_inverse_status(certificate),
        copied_polynomial_column =
            _case008_d14_ltp_copied_polynomial_column_status(certificate),
        modified_witness_coefficient =
            _case008_d14_ltp_modified_witness_status(certificate),
        copied_witness_from_different_column =
            _case008_d14_ltp_copied_witness_status(certificate),
    )
    controls_pass = controls.nonunimodular_source == :rejected &&
        all(status -> status != :ok, values(controls)[2:end])
    return controls, controls_pass
end

function _case008_d14_ltp_build_certificate(fixture)
    fixture_validation = ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture)
    fixture_validation == :ok ||
        throw(ArgumentError("invalid case_008 d14 fixture: $(fixture_validation)"))
    selected_generator = _case008_d14_ltp_selected_generator(fixture.ring)
    noether = Suslin._laurent_noether_certificate(
        fixture.failing_column,
        CASE008_D14_LTP_SELECTED_ENTRY,
        selected_generator,
    )
    certificate = Suslin._laurent_to_polynomial_certificate(
        fixture.failing_column,
        noether,
        CASE008_D14_LTP_SELECTED_ENTRY,
        selected_generator,
    )
    return certificate
end

function case008_d14_laurent_to_polynomial_conversion_report(;
    timeout_seconds = 300,
)
    return case008_d14_laurent_to_polynomial_conversion_report(
        _case008_d14_ltp_default_fixture();
        timeout_seconds,
    )
end

function case008_d14_laurent_to_polynomial_conversion_report(
    fixture;
    timeout_seconds = 300,
)
    haskey(CASE008_D14_LTP_REPORT_CACHE, fixture) &&
        return CASE008_D14_LTP_REPORT_CACHE[fixture].report
    timeout_seconds > 0 ||
        throw(ArgumentError("timeout_seconds must be positive"))

    start_time = time()
    certificate = _case008_d14_ltp_build_certificate(fixture)
    elapsed_seconds = time() - start_time

    conversion_status = Suslin._validate_laurent_to_polynomial_certificate(certificate)
    replay = certificate.replay
    negative_controls, negative_controls_pass =
        _case008_d14_ltp_negative_controls(fixture, certificate)
    report = (;
        case_id = fixture.case_id,
        source_dimension = length(fixture.failing_column),
        ring_generators = _case008_d14_ltp_ring_generator_names(fixture.ring),
        source_boundary = :ToricBuilderCase008D14ColumnBoundary,
        selected_entry = CASE008_D14_LTP_SELECTED_ENTRY,
        selected_generator = CASE008_D14_LTP_SELECTED_GENERATOR_NAME,
        conversion_status,
        conversion_replay = replay.overall_ok ? :pass : :fail,
        conversion_factors = length(certificate.conversion_factors),
        source_terms =
            Suslin._laurent_to_polynomial_column_term_count(fixture.failing_column),
        polynomial_terms =
            Suslin._laurent_to_polynomial_column_term_count(certificate.polynomial_column),
        polynomialization_power = certificate.polynomialization_power,
        polynomial_unimodular = replay.polynomial_unimodular,
        polynomial_witness_replay =
            replay.polynomial_witness_replay_ok ? :pass : :fail,
        polynomial_witness_terms =
            certificate.polynomial_witness_metadata.witness_terms,
        polynomial_witness_provenance =
            certificate.polynomial_witness_metadata.provenance,
        negative_controls,
        negative_controls_pass,
        elapsed_seconds,
        timeout_seconds,
        final_status =
            (conversion_status == :ok && negative_controls_pass && elapsed_seconds <= timeout_seconds) ?
                :ok : :timed_out_or_failed,
        status = CASE008_D14_LTP_REPORT_STATUS,
    )
    CASE008_D14_LTP_REPORT_CACHE[fixture] = (; fixture, report)
    return report
end

function _case008_d14_ltp_report_has_required_fields(report)::Bool
    return all(field -> hasproperty(report, field), CASE008_D14_LTP_REPORT_FIELDS)
end

function validate_case008_d14_laurent_to_polynomial_conversion_report(report)::Symbol
    return validate_case008_d14_laurent_to_polynomial_conversion_report(
        report,
        _case008_d14_ltp_default_fixture(),
    )
end

function validate_case008_d14_laurent_to_polynomial_conversion_report(
    report,
    fixture,
)::Symbol
    try
        ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok ||
            return :invalid_fixture
        _case008_d14_ltp_report_has_required_fields(report) ||
            return :missing_report_fields
        report.status == CASE008_D14_LTP_REPORT_STATUS || return :wrong_status
        report.case_id == "case_008" || return :wrong_case
        report.source_dimension == 14 || return :wrong_source_dimension
        report.ring_generators == ("u", "v") || return :wrong_ring_generators
        report.source_boundary == :ToricBuilderCase008D14ColumnBoundary ||
            return :wrong_source_boundary
        report.selected_entry == CASE008_D14_LTP_SELECTED_ENTRY ||
            return :wrong_selected_entry
        report.selected_generator == CASE008_D14_LTP_SELECTED_GENERATOR_NAME ||
            return :wrong_selected_generator
        report.conversion_status == :ok || return :conversion_failed
        report.conversion_replay == :pass || return :conversion_replay_failed
        report.conversion_factors == CASE008_D14_LTP_EXPECTED_CONVERSION_FACTORS ||
            return :wrong_conversion_factor_count
        report.source_terms == CASE008_D14_LTP_EXPECTED_SOURCE_TERMS ||
            return :wrong_source_term_count
        report.polynomial_terms == CASE008_D14_LTP_EXPECTED_POLYNOMIAL_TERMS ||
            return :wrong_polynomial_term_count
        report.polynomialization_power ==
            CASE008_D14_LTP_EXPECTED_POLYNOMIALIZATION_POWER ||
            return :wrong_polynomialization_power
        report.polynomial_unimodular === true || return :not_polynomial_unimodular
        report.polynomial_witness_replay == :pass ||
            return :polynomial_witness_replay_failed
        report.polynomial_witness_terms isa Integer ||
            return :wrong_polynomial_witness_terms
        report.polynomial_witness_terms >= 0 ||
            return :wrong_polynomial_witness_terms
        report.polynomial_witness_provenance isa Symbol ||
            return :wrong_polynomial_witness_provenance
        report.negative_controls_pass === true || return :negative_controls_failed
        report.elapsed_seconds isa Real || return :wrong_elapsed_seconds
        isfinite(report.elapsed_seconds) || return :wrong_elapsed_seconds
        report.elapsed_seconds >= 0 || return :wrong_elapsed_seconds
        report.timeout_seconds isa Real || return :wrong_timeout_seconds
        isfinite(report.timeout_seconds) || return :wrong_timeout_seconds
        report.elapsed_seconds <= report.timeout_seconds || return :timed_out
        report.final_status == :ok || return :wrong_final_status
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_report
    end
end

function _case008_d14_ltp_value_text(value)
    value isa Symbol && return String(value)
    value isa AbstractFloat && return string(round(value; digits = 3))
    value isa Tuple && return join(_case008_d14_ltp_value_text.(value), ", ")
    return string(value)
end

function render_case008_d14_laurent_to_polynomial_conversion_markdown(report)
    io = IOBuffer()
    println(io, "# case_008 d14 LaurentToPoly Conversion Report")
    println(io)
    for field in CASE008_D14_LTP_REPORT_FIELDS
        field == :negative_controls && continue
        println(io, "$(field) = $(_case008_d14_ltp_value_text(getproperty(report, field)))")
    end
    println(io)
    println(io, "## Negative Controls")
    println(io)
    println(io, "| control | status |")
    println(io, "| --- | --- |")
    for (name, status) in pairs(report.negative_controls)
        println(io, "| $(name) | $(_case008_d14_ltp_value_text(status)) |")
    end
    return String(take!(io))
end

function write_case008_d14_laurent_to_polynomial_conversion_report(
    path::AbstractString;
    report = case008_d14_laurent_to_polynomial_conversion_report(),
)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, render_case008_d14_laurent_to_polynomial_conversion_markdown(report))
    end
    return path
end

function _case008_d14_ltp_parse_timeout_seconds(raw::AbstractString)
    value = try
        parse(Float64, raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("--timeout-seconds must be a positive number, got $(repr(raw))"))
    end
    value > 0 || throw(ArgumentError("--timeout-seconds must be positive, got $(raw)"))
    return value
end

function _case008_d14_ltp_parse_args(args)
    output = joinpath(tempdir(), "case008-d14-ltp-conversion.md")
    timeout_seconds = 300.0
    for arg in args
        if startswith(arg, "--output=")
            output = arg[length("--output=")+1:end]
        elseif startswith(arg, "--timeout-seconds=")
            timeout_seconds =
                _case008_d14_ltp_parse_timeout_seconds(arg[length("--timeout-seconds=")+1:end])
        else
            throw(ArgumentError("unsupported argument: $(arg)"))
        end
    end
    return (; output, timeout_seconds)
end

function main(args = ARGS)
    options = _case008_d14_ltp_parse_args(args)
    report = case008_d14_laurent_to_polynomial_conversion_report(;
        timeout_seconds = options.timeout_seconds,
    )
    path = write_case008_d14_laurent_to_polynomial_conversion_report(
        options.output;
        report,
    )
    validation = validate_case008_d14_laurent_to_polynomial_conversion_report(report)
    validation == :ok ||
        throw(ArgumentError("case_008 d14 LaurentToPoly conversion report validation failed: $(validation)"))
    println("wrote case_008 d14 LaurentToPoly conversion report to $(path)")
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
