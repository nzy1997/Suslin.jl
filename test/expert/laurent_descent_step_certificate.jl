using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_elementary_move_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_elementary_move_search.jl"))
end

const LAURENT_DESCENT_STEP_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :operation,
    :before_profile,
    :after_profile,
    :before_measure,
    :after_measure,
    :status,
    :replay_status,
    :measure_relation,
)

const LAURENT_DESCENT_STEP_OPERATION_FIELDS = (
    :family,
    :target_index,
    :source_index,
    :coefficient,
    :exponent,
    :ring_generators,
)

function laurent_descent_step_profile(column, R; case_id)
    generator_names = _ring_generator_names(R)
    support_summary = _newton_support_summary(column, generator_names)
    term_counts = support_summary.per_entry_support_counts
    return (;
        case_id,
        dimension = length(column),
        ring_generators = generator_names,
        nonzero_entries = count(!iszero, column),
        max_entry_terms = maximum(term_counts; init = 0),
        entry_term_counts = term_counts,
        valuation_ranges = _valuation_ranges(
            support_summary.whole_column_bounds,
            generator_names,
        ),
        newton_support_summary = support_summary,
        leading_monomial_candidates = _leading_monomial_candidates(
            column,
            generator_names,
        ),
        candidate_measure_families = CASE008_D14_PROFILE_MEASURE_FAMILIES,
        status = :profile_only,
    )
end

function _laurent_descent_step_operation_status(operation, n::Int, R)::Symbol
    _has_required_fields(operation, LAURENT_DESCENT_STEP_OPERATION_FIELDS) ||
        return :malformed_operation
    operation.family == :entry_addition || return :malformed_operation
    operation.ring_generators == _ring_generator_names(R) ||
        return :wrong_ring_generators
    try
        target = _checked_entry_index(operation.target_index, n, "target_index")
        source = _checked_entry_index(operation.source_index, n, "source_index")
        target != source || return :malformed_operation
        _checked_exponent_pair(operation.exponent)
        R(operation.coefficient)
    catch err
        err isa InterruptException && rethrow()
        return :malformed_operation
    end
    return :ok
end

function _laurent_descent_step_certificate_from_replay(
    column,
    R,
    before_profile,
    operation;
    require_strict::Bool = true,
)
    hasproperty(before_profile, :case_id) ||
        throw(ArgumentError("before_profile must include case_id"))
    operation_status =
        _laurent_descent_step_operation_status(operation, length(column), R)
    operation_status == :ok ||
        throw(ArgumentError("invalid Laurent descent operation: $(operation_status)"))

    expected_before_profile = laurent_descent_step_profile(
        column,
        R;
        case_id = before_profile.case_id,
    )
    before_profile == expected_before_profile ||
        throw(ArgumentError("before_profile is stale for the input column"))
    before_measure = _case008_d14_measure_from_column(
        column,
        R;
        case_id = before_profile.case_id,
    )
    after_column = replay_laurent_elementary_entry_addition(
        column,
        R,
        operation,
    )
    after_profile = laurent_descent_step_profile(
        after_column,
        R;
        case_id = before_profile.case_id,
    )
    after_measure = _case008_d14_measure_from_column(
        after_column,
        R;
        case_id = before_profile.case_id,
    )
    relation = strictly_decreases_laurent_measure(
        before_measure,
        after_measure,
    ) ? :strict_decrease : :not_strict_decrease
    require_strict && relation == :not_strict_decrease &&
        throw(ArgumentError("operation does not strictly decrease the measure"))
    return (;
        case_id = before_profile.case_id,
        dimension = length(column),
        ring_generators = _ring_generator_names(R),
        operation,
        before_profile,
        after_profile,
        before_measure,
        after_measure,
        status = :descent_step_certificate,
        replay_status = :ok,
        measure_relation = relation,
    )
end

function laurent_descent_step_certificate(column, R, before_profile, operation)
    return _laurent_descent_step_certificate_from_replay(
        column,
        R,
        before_profile,
        operation;
        require_strict = true,
    )
end

const CASE008_D14_RECORDED_DESCENT_OPERATION = (;
    family = :entry_addition,
    target_index = 1,
    source_index = 2,
    coefficient = 1,
    exponent = (-1, 1),
    ring_generators = ("u", "v"),
)

function case008_d14_laurent_descent_step_certificate(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    profile = case008_d14_laurent_descent_profile(fixture)
    return laurent_descent_step_certificate(
        fixture.failing_column,
        fixture.ring,
        profile,
        CASE008_D14_RECORDED_DESCENT_OPERATION,
    )
end

function validate_laurent_descent_step_certificate(cert, column, R)::Symbol
    try
        _has_required_fields(cert, LAURENT_DESCENT_STEP_CERTIFICATE_FIELDS) ||
            return :missing_certificate_fields
        cert.status == :descent_step_certificate || return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.measure_relation == :strict_decrease ||
            return :wrong_measure_relation
        cert.dimension == length(column) || return :wrong_dimension

        ring_generators = _ring_generator_names(R)
        cert.ring_generators == ring_generators ||
            return :wrong_ring_generators

        operation_status = _laurent_descent_step_operation_status(
            cert.operation,
            length(column),
            R,
        )
        operation_status == :ok || return operation_status

        hasproperty(cert.before_profile, :case_id) ||
            return :missing_profile_fields
        cert.case_id == cert.before_profile.case_id || return :wrong_case

        expected_before_profile = laurent_descent_step_profile(
            column,
            R;
            case_id = cert.case_id,
        )
        cert.before_profile == expected_before_profile ||
            return :stale_before_profile

        expected_before_measure = _case008_d14_measure_from_column(
            column,
            R;
            case_id = cert.case_id,
        )
        cert.before_measure == expected_before_measure ||
            return :stale_before_measure

        after_column = replay_laurent_elementary_entry_addition(
            column,
            R,
            cert.operation,
        )
        expected_after_profile = laurent_descent_step_profile(
            after_column,
            R;
            case_id = cert.case_id,
        )
        cert.after_profile == expected_after_profile ||
            return :stale_after_profile

        expected_after_measure = _case008_d14_measure_from_column(
            after_column,
            R;
            case_id = cert.case_id,
        )
        cert.after_measure == expected_after_measure ||
            return :stale_after_measure

        strictly_decreases_laurent_measure(
            expected_before_measure,
            expected_after_measure,
        ) || return :not_strict_decrease
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :operation_replay_failed
    end
end

@testset "Laurent descent-step certificate shell" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/laurent_descent_step_certificate.jl\"",
        runtests,
    )

    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [u + v, u + v + one(R)]
    before_profile = laurent_descent_step_profile(
        column,
        R;
        case_id = "synthetic",
    )
    operation = (;
        family = :entry_addition,
        target_index = 2,
        source_index = 1,
        coefficient = 1,
        exponent = (0, 0),
        ring_generators = ("u", "v"),
    )

    cert = laurent_descent_step_certificate(
        column,
        R,
        before_profile,
        operation,
    )

    @test cert.case_id == "synthetic"
    @test cert.dimension == 2
    @test cert.ring_generators == ("u", "v")
    @test cert.operation == operation
    @test cert.status == :descent_step_certificate
    @test cert.replay_status == :ok
    @test cert.measure_relation == :strict_decrease
    @test strictly_decreases_laurent_measure(
        cert.before_measure,
        cert.after_measure,
    )
    @test validate_laurent_descent_step_certificate(cert, column, R) == :ok

    tampered_operation = merge(operation, (; exponent = (1, 0)))
    tampered_cert = merge(cert, (; operation = tampered_operation))
    @test validate_laurent_descent_step_certificate(
        tampered_cert,
        column,
        R,
    ) == :stale_after_profile

    stale_after_profile = merge(
        cert.after_profile,
        (; nonzero_entries = cert.after_profile.nonzero_entries + 1),
    )
    stale_after_cert = merge(cert, (; after_profile = stale_after_profile))
    @test validate_laurent_descent_step_certificate(
        stale_after_cert,
        column,
        R,
    ) == :stale_after_profile

    equal_operation = merge(operation, (; coefficient = 0))
    equal_cert = _laurent_descent_step_certificate_from_replay(
        column,
        R,
        before_profile,
        equal_operation;
        require_strict = false,
    )
    equal_claim = merge(equal_cert, (; measure_relation = :strict_decrease))
    @test validate_laurent_descent_step_certificate(
        equal_claim,
        column,
        R,
    ) == :not_strict_decrease

    wrong_ring_cert = merge(cert, (; ring_generators = ("v", "u")))
    @test validate_laurent_descent_step_certificate(
        wrong_ring_cert,
        column,
        R,
    ) == :wrong_ring_generators

    malformed_operation = (;
        family = :entry_addition,
        target_index = 2,
        source_index = 1,
        coefficient = 1,
        ring_generators = ("u", "v"),
    )
    malformed_cert = merge(cert, (; operation = malformed_operation))
    @test validate_laurent_descent_step_certificate(
        malformed_cert,
        column,
        R,
    ) == :malformed_operation

    wrong_status = merge(cert, (; status = :profile_only))
    @test validate_laurent_descent_step_certificate(
        wrong_status,
        column,
        R,
    ) == :wrong_status
end

@testset "case_008 d=14 Laurent descent-step certificate" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    cert = case008_d14_laurent_descent_step_certificate(fixture)
    @test cert.case_id == "case_008"
    @test cert.dimension == 14
    @test cert.ring_generators == ("u", "v")
    @test cert.operation.family == :entry_addition
    @test cert.operation.target_index == 1
    @test cert.operation.source_index == 2
    @test cert.operation.coefficient == 1
    @test cert.operation.exponent == (-1, 1)
    @test cert.operation.ring_generators == ("u", "v")
    @test cert.status == :descent_step_certificate
    @test cert.replay_status == :ok
    @test cert.measure_relation == :strict_decrease
    @test validate_laurent_descent_step_certificate(
        cert,
        fixture.failing_column,
        fixture.ring,
    ) == :ok
end
