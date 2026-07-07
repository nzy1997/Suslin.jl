using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl"))

const INTERNAL_D14_OPERATION = (;
    family = :entry_addition,
    target_index = 1,
    source_index = 2,
    coefficient = 1,
    exponent = (-1, 1),
    ring_generators = ("u", "v"),
)

function _internal_descent_certificate(column, R, operation; case_id = "case_008")
    before = Suslin._laurent_descent_measure_from_column(column, R; case_id)
    after_column = Suslin._replay_laurent_elementary_entry_addition(column, R, operation)
    after = Suslin._laurent_descent_measure_from_column(after_column, R; case_id)
    return (;
        case_id,
        dimension = length(column),
        ring_generators = Tuple(string.(gens(R))),
        operation,
        before_measure = before,
        after_measure = after,
        status = :descent_step_certificate,
        replay_status = :ok,
        measure_relation = Suslin._strictly_decreases_laurent_measure(before, after) ?
            :strict_decrease : :not_strict_decrease,
    )
end

@testset "internal Laurent descent measure helpers" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin("\"internal/laurent_descent_measure_helpers.jl\"", runtests)

    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    measure = Suslin._laurent_descent_measure_from_column(
        fixture.failing_column,
        fixture.ring;
        case_id = fixture.case_id,
    )
    @test measure.case_id == "case_008"
    @test measure.dimension == 14
    @test measure.ring_generators == ("u", "v")
    @test measure.status == :measure_contract
    @test measure.order == :lexicographic_minimize
    @test measure.components == (
        :whole_support_count,
        :max_entry_terms,
        :valuation_span,
        :leading_exponent,
        :leading_entry_index,
    )
    @test measure.whole_support_count == 7387
    @test measure.max_entry_terms == 3734
    @test measure.valuation_span == (97, 93)
    @test measure.leading_exponent == (49, -5)
    @test measure.leading_entry_index == 10

    after_column = Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        INTERNAL_D14_OPERATION,
    )
    after_measure = Suslin._laurent_descent_measure_from_column(
        after_column,
        fixture.ring;
        case_id = fixture.case_id,
    )
    @test Suslin._strictly_decreases_laurent_measure(measure, after_measure)

    cert = _internal_descent_certificate(
        fixture.failing_column,
        fixture.ring,
        INTERNAL_D14_OPERATION;
        case_id = fixture.case_id,
    )
    @test Suslin._validate_laurent_descent_step_certificate(
        cert,
        fixture.failing_column,
        fixture.ring,
    ) == :ok

    swapped_cert = merge(cert, (; ring_generators = ("v", "u")))
    @test Suslin._validate_laurent_descent_step_certificate(
        swapped_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :wrong_ring_generators

    swapped_operation = merge(INTERNAL_D14_OPERATION, (; ring_generators = ("v", "u")))
    swapped_operation_cert = merge(cert, (; operation = swapped_operation))
    @test Suslin._validate_laurent_descent_step_certificate(
        swapped_operation_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :wrong_ring_generators

    missing_ring_generators_operation = Base.structdiff(
        INTERNAL_D14_OPERATION,
        (; ring_generators = INTERNAL_D14_OPERATION.ring_generators),
    )
    missing_ring_generators_cert = merge(
        cert,
        (; operation = missing_ring_generators_operation),
    )
    @test Suslin._validate_laurent_descent_step_certificate(
        missing_ring_generators_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :malformed_operation

    bad_target = merge(INTERNAL_D14_OPERATION, (; target_index = 0))
    @test_throws ArgumentError Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        bad_target,
    )
    bad_target_cert = merge(cert, (; operation = bad_target))
    @test Suslin._validate_laurent_descent_step_certificate(
        bad_target_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :malformed_operation

    bad_source = merge(INTERNAL_D14_OPERATION, (; source_index = 0))
    @test_throws ArgumentError Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        bad_source,
    )
    bad_source_cert = merge(cert, (; operation = bad_source))
    @test Suslin._validate_laurent_descent_step_certificate(
        bad_source_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :malformed_operation

    equal_indices = merge(INTERNAL_D14_OPERATION, (; source_index = 1))
    @test_throws ArgumentError Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        equal_indices,
    )
    equal_indices_cert = merge(cert, (; operation = equal_indices))
    @test Suslin._validate_laurent_descent_step_certificate(
        equal_indices_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :malformed_operation

    stale_before = merge(
        cert,
        (;
            before_measure = merge(
                cert.before_measure,
                (; whole_support_count = cert.before_measure.whole_support_count - 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_descent_step_certificate(
        stale_before,
        fixture.failing_column,
        fixture.ring,
    ) == :stale_before_measure

    zero_operation = merge(INTERNAL_D14_OPERATION, (; coefficient = 0))
    zero_cert = _internal_descent_certificate(
        fixture.failing_column,
        fixture.ring,
        zero_operation;
        case_id = fixture.case_id,
    )
    zero_claim = merge(zero_cert, (; measure_relation = :strict_decrease))
    @test Suslin._validate_laurent_descent_step_certificate(
        zero_claim,
        fixture.failing_column,
        fixture.ring,
    ) == :not_strict_decrease

    stale_after_operation = merge(INTERNAL_D14_OPERATION, (; exponent = (0, 0)))
    stale_after_cert = merge(cert, (; operation = stale_after_operation))
    @test Suslin._validate_laurent_descent_step_certificate(
        stale_after_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :stale_after_measure

    bad_replay_column = Any[fixture.failing_column...]
    bad_replay_column[1] = nothing
    @test Suslin._validate_laurent_descent_step_certificate(
        cert,
        bad_replay_column,
        fixture.ring,
    ) == :operation_replay_failed

    P, _ = polynomial_ring(QQ, ["u", "v"])
    @test Suslin._laurent_descent_step_diagnostic_certificate(
        [one(P) for _ in 1:length(fixture.failing_column)],
        P,
    ) === nothing
end
