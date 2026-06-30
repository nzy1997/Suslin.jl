using Test
using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")

function _context_as_namedtuple(context)
    names = propertynames(context)
    return NamedTuple{names}(Tuple(getproperty(context, name) for name in names))
end

function _corrupt_context(context, updates)
    return Suslin.SL3RealizationInputContext(
        values(merge(_context_as_namedtuple(context), updates))...,
    )
end

function _corruption_is_rejected(context, updates)
    try
        corrupted_context = _corrupt_context(context, updates)
        return !Suslin._verify_sl3_realization_input_context(corrupted_context)
    catch err
        return err isa ArgumentError
    end
end

function _diagnostic_mentions(diagnostic, field::Symbol, target)
    hasproperty(diagnostic, field) || return false
    value = getproperty(diagnostic, field)
    return value === target || (value isa AbstractArray || value isa Tuple || value isa Set ?
        target in value :
        false)
end

function _catalog_metadata(entry)
    return (;
        fixture_id = entry.id,
        role = entry.role,
        expected_status = entry.expected_status,
        consumer_issue_ids = entry.consumer_issue_ids,
    )
end

if !isdefined(Main, :ParkWoodburnSL3DriverFixtureCatalog)
    include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
end

@testset "Park-Woodburn SL3 driver context contract" begin
    catalog = ParkWoodburnSL3DriverFixtureCatalog.catalog()
    entries = ParkWoodburnSL3DriverFixtureCatalog.cases_by_id()
    negative = Dict(entry.id => entry for entry in catalog.negative_controls)

    fast_entry = entries["sl3-driver-univariate-fast-local-qq"]
    fast_context = Suslin._sl3_realization_input_context(
        fast_entry.matrix;
        selected_variable = fast_entry.selected_variable,
        catalog_metadata = _catalog_metadata(fast_entry),
        local_form_witness = fast_entry.local_form_witness,
    )
    @test fast_context.support_status == :supported
    @test fast_context.evidence_status == :replayable
    @test fast_context.local_form_status == :replayed
    @test fast_context.determinant_status == :one
    @test fast_context.exact_field_status == :supported
    @test Suslin._verify_sl3_realization_input_context(fast_context)

    legacy_entry = entries["sl3-driver-legacy-quillen-patched-substitution-qq"]
    legacy_context = Suslin._sl3_realization_input_context(
        legacy_entry.matrix;
        selected_variable = legacy_entry.selected_variable,
        catalog_metadata = _catalog_metadata(legacy_entry),
        quillen_murthy_metadata = legacy_entry.upstream_evidence,
    )
    @test legacy_context.support_status == :staged
    @test legacy_context.evidence_status == :partial
    @test legacy_context.quillen_murthy_status == :recorded
    @test _diagnostic_mentions(legacy_context.staged_diagnostic, :partial_evidence, :quillen_murthy) ||
          _diagnostic_mentions(legacy_context.staged_diagnostic, :missing_evidence, :quillen_murthy)
    @test Suslin._verify_sl3_realization_input_context(legacy_context)

    multivariate_entry = entries["sl3-driver-multivariate-monic-special-form-qq"]
    multivariate_context = Suslin._sl3_realization_input_context(
        multivariate_entry.matrix;
        selected_variable = multivariate_entry.selected_variable,
        catalog_metadata = _catalog_metadata(multivariate_entry),
        local_form_witness = multivariate_entry.local_form_witness,
    )
    @test multivariate_context.support_status == :supported
    @test multivariate_context.evidence_status == :replayable
    @test multivariate_context.local_form_status == :replayed
    @test Suslin._verify_sl3_realization_input_context(multivariate_context)

    staged_entry = entries["sl3-driver-det-one-no-witness-staged-qq"]
    staged_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = _catalog_metadata(staged_entry),
    )
    @test staged_context.support_status == :staged
    @test staged_context.evidence_status == :missing
    @test Suslin._verify_sl3_realization_input_context(staged_context)

    fake_quillen_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = merge(_catalog_metadata(staged_entry), (; expected_status = :staged,)),
        quillen_murthy_metadata = (; fixture_id = :fake, matrix = staged_entry.matrix),
    )
    @test fake_quillen_context.support_status == :staged
    @test fake_quillen_context.evidence_status == :partial
    @test fake_quillen_context.quillen_murthy_status == :recorded
    @test Suslin._verify_sl3_realization_input_context(fake_quillen_context)

    variable_change_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "synthetic-variable-change-replay"),
        variable_change_metadata = (;
            replay_id = "variable-change-replay",
            target_matrix = staged_entry.matrix,
            replay_steps = ((; kind = :supplied_variable_change_replay),),
        ),
    )
    @test variable_change_context.support_status == :supported
    @test variable_change_context.evidence_status == :replayable
    @test variable_change_context.variable_change_status == :replayed
    @test Suslin._verify_sl3_realization_input_context(variable_change_context)

    variable_change_recorded_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = _catalog_metadata(staged_entry),
        variable_change_metadata = (;
            replay_id = "variable-change-recorded",
            target_matrix = staged_entry.matrix,
        ),
    )
    @test variable_change_recorded_context.support_status == :staged
    @test variable_change_recorded_context.evidence_status == :partial
    @test variable_change_recorded_context.variable_change_status == :recorded
    @test Suslin._verify_sl3_realization_input_context(variable_change_recorded_context)

    normality_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "synthetic-normality-replay"),
        normality_conjugation_metadata = (;
            replay_id = "normality-conjugation-replay",
            target_matrix = staged_entry.matrix,
            replay_steps = ((; kind = :supplied_normality_conjugation_replay),),
        ),
    )
    @test normality_context.support_status == :supported
    @test normality_context.evidence_status == :replayable
    @test normality_context.normality_conjugation_status == :replayed
    @test Suslin._verify_sl3_realization_input_context(normality_context)

    normality_recorded_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = _catalog_metadata(staged_entry),
        normality_conjugation_metadata = (;
            replay_id = "normality-conjugation-recorded",
            target_matrix = staged_entry.matrix,
        ),
    )
    @test normality_recorded_context.support_status == :staged
    @test normality_recorded_context.evidence_status == :partial
    @test normality_recorded_context.normality_conjugation_status == :recorded
    @test Suslin._verify_sl3_realization_input_context(normality_recorded_context)

    @test_throws ArgumentError Suslin._sl3_realization_input_context(
        negative["sl3-driver-negative-det-not-one"].matrix;
        selected_variable =
            negative["sl3-driver-negative-det-not-one"].selected_variable.generator,
        catalog_metadata =
            _catalog_metadata(negative["sl3-driver-negative-det-not-one"]),
    )
    @test_throws ArgumentError Suslin._sl3_realization_input_context(
        negative["sl3-driver-negative-unsupported-coefficient-ring"].matrix;
        selected_variable =
            negative["sl3-driver-negative-unsupported-coefficient-ring"].selected_variable.generator,
        catalog_metadata =
            _catalog_metadata(negative["sl3-driver-negative-unsupported-coefficient-ring"]),
        local_form_witness =
            negative["sl3-driver-negative-unsupported-coefficient-ring"].local_form_witness,
    )
    @test_throws ArgumentError Suslin._sl3_realization_input_context(
        negative["sl3-driver-negative-selected-variable-not-generator"].matrix;
        selected_variable =
            negative["sl3-driver-negative-selected-variable-not-generator"].selected_variable.generator,
        catalog_metadata =
            _catalog_metadata(negative["sl3-driver-negative-selected-variable-not-generator"]),
    )
    @test_throws ArgumentError Suslin._sl3_realization_input_context(
        fast_entry.matrix;
        selected_variable = merge(fast_entry.selected_variable, (; index = 2,)),
        catalog_metadata = _catalog_metadata(fast_entry),
        local_form_witness = fast_entry.local_form_witness,
    )
    @test_throws ArgumentError Suslin._sl3_realization_input_context(
        negative["sl3-driver-negative-supported-without-witness"].matrix;
        selected_variable =
            negative["sl3-driver-negative-supported-without-witness"].selected_variable.generator,
        catalog_metadata =
            _catalog_metadata(negative["sl3-driver-negative-supported-without-witness"]),
    )

    @test _corruption_is_rejected(
        fast_context,
        (; selected_variable = one(parent(fast_entry.selected_variable.generator))),
    )
    @test _corruption_is_rejected(fast_context, (; determinant_status = :not_one))
    @test _corruption_is_rejected(fast_context, (; ring_profile = :unsupported))
    @test _corruption_is_rejected(fast_context, (; evidence_status = :missing))
    @test _corruption_is_rejected(
        legacy_context,
        (; staged_diagnostic = merge(
            legacy_context.staged_diagnostic,
            (; partial_evidence = ("corrupted_evidence",)),
        )),
    )
end
