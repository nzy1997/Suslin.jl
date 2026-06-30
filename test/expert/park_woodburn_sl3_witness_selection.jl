using Test
using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")

function _selection_as_namedtuple(selection)
    names = propertynames(selection)
    return NamedTuple{names}(Tuple(getproperty(selection, name) for name in names))
end

function _corrupt_selection(selection, updates)
    return Suslin.SL3LocalFormWitnessSelection(
        values(merge(_selection_as_namedtuple(selection), updates))...,
    )
end

function _selection_corruption_is_rejected(selection, updates)
    try
        corrupted = _corrupt_selection(selection, updates)
        return !Suslin._verify_sl3_local_form_witness_selection(corrupted)
    catch err
        return err isa ArgumentError
    end
end

function _catalog_metadata(entry; expected_status = entry.expected_status)
    return (;
        fixture_id = entry.id,
        role = entry.role,
        expected_status,
        consumer_issue_ids = entry.consumer_issue_ids,
    )
end

function _special_form_matrix(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

if !isdefined(Main, :ParkWoodburnSL3DriverFixtureCatalog)
    include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
end

@testset "Park-Woodburn SL3 local witness selection" begin
    entries = ParkWoodburnSL3DriverFixtureCatalog.cases_by_id()

    fast_entry = entries["sl3-driver-univariate-fast-local-qq"]
    fast_context = Suslin._sl3_realization_input_context(
        fast_entry.matrix;
        selected_variable = fast_entry.selected_variable,
        catalog_metadata = _catalog_metadata(fast_entry),
    )
    fast_selection = Suslin._select_sl3_local_form_witness(fast_context)
    @test fast_selection.support_status == :supported
    @test fast_selection.replay_status == :replayed
    @test fast_selection.witness_source == :already_special_form
    @test fast_selection.selected_variable == fast_entry.selected_variable.generator
    @test fast_selection.selected_variable_index == fast_entry.selected_variable.index
    @test fast_selection.entries ==
          (; p = fast_entry.matrix[1, 1], q = fast_entry.matrix[1, 2],
             r = fast_entry.matrix[2, 1], s = fast_entry.matrix[2, 2])
    @test fast_selection.monicity_witness.is_monic == true
    @test fast_selection.monicity_witness.variable == fast_entry.selected_variable.generator
    @test Suslin._verify_sl3_realization_input_context(fast_context)
    @test Suslin._verify_sl3_local_form_witness_selection(fast_selection)

    staged_entry = entries["sl3-driver-det-one-no-witness-staged-qq"]
    R = base_ring(staged_entry.matrix)
    X, u, v = gens(R)
    p = X + u * v + one(R)
    q = one(R)
    r = X + u * v
    s = one(R)
    multivariate_target = _special_form_matrix(R, p, q, r, s)
    @test det(multivariate_target) == one(R)
    multivariate_context = Suslin._sl3_realization_input_context(
        multivariate_target;
        selected_variable = (; name = "X", generator = X, index = 1, status = :passes),
        catalog_metadata = (; fixture_id = "issue-236-multivariate-special-form"),
        local_form_witness = (; entries = (; p, q, r, s)),
    )
    multivariate_selection = Suslin._select_sl3_local_form_witness(multivariate_context)
    @test multivariate_selection.support_status == :supported
    @test multivariate_selection.entries == (; p, q, r, s)
    @test multivariate_selection.local_form_matrix == multivariate_target
    @test multivariate_selection.monicity_witness.degree == degree(p, 1)
    @test multivariate_selection.monicity_witness.leading_coefficient == one(R)
    @test Suslin._verify_sl3_realization_input_context(multivariate_context)
    @test Suslin._verify_sl3_local_form_witness_selection(multivariate_selection)

    staged_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = _catalog_metadata(staged_entry),
    )
    staged_selection = Suslin._select_sl3_local_form_witness(staged_context)
    @test staged_selection.support_status == :staged
    @test staged_selection.replay_status == :missing
    @test staged_selection.entries === nothing
    @test occursin(
        "missing supported local-form witness",
        staged_selection.staged_diagnostic.reason,
    )
    @test Suslin._verify_sl3_realization_input_context(staged_context)
    @test Suslin._verify_sl3_local_form_witness_selection(staged_selection)

    source_matrix = staged_entry.matrix
    variable_change_metadata = (;
        replay_id = "issue-236-variable-change-replay",
        source_matrix,
        selected_variable = staged_entry.selected_variable.generator,
        local_form_matrix = multivariate_target,
        replay_steps = ((; kind = :supplied_variable_change, name = :identity_catalog_replay),),
    )
    variable_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-supplied-variable-change"),
        variable_change_metadata,
    )
    variable_selection = Suslin._select_sl3_local_form_witness(variable_context)
    @test variable_selection.support_status == :supported
    @test variable_selection.replay_status == :replayed
    @test variable_selection.witness_source == :variable_change
    @test variable_selection.variable_change_status == :replayed
    @test variable_selection.entries == (; p, q, r, s)
    @test variable_selection.variable_change_metadata == variable_change_metadata
    @test Suslin._verify_sl3_realization_input_context(variable_context)
    @test Suslin._verify_sl3_local_form_witness_selection(variable_selection)

    conflicting_p = X + u * v + R(2)
    conflicting_r = conflicting_p - one(R)
    conflicting_target = _special_form_matrix(R, conflicting_p, q, conflicting_r, s)
    @test det(conflicting_target) == one(R)
    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        multivariate_context;
        variable_change_metadata = (;
            replay_id = "issue-236-conflicting-variable-change-replay",
            source_matrix = multivariate_target,
            selected_variable = X,
            local_form_matrix = conflicting_target,
            replay_steps = ((; kind = :conflicting_variable_change),),
        ),
    )

    normality_shell_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-normality-shell"),
        normality_conjugation_metadata = (;
            replay_id = "normality-shell-without-local-form",
            source_matrix,
            replay_steps = ((; kind = :normality_shell),),
        ),
    )
    normality_shell_selection =
        Suslin._select_sl3_local_form_witness(normality_shell_context)
    @test normality_shell_selection.support_status == :staged
    @test normality_shell_selection.normality_conjugation_status == :recorded
    @test :normality_conjugation in
          normality_shell_selection.staged_diagnostic.partial_evidence
    @test Suslin._verify_sl3_realization_input_context(normality_shell_context)
    @test Suslin._verify_sl3_local_form_witness_selection(normality_shell_selection)

    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        multivariate_context;
        selected_variable = X + one(R),
    )

    nonmonic_p = 2 * X + u * v + one(R)
    nonmonic_r = nonmonic_p - one(R)
    nonmonic_target = _special_form_matrix(R, nonmonic_p, q, nonmonic_r, s)
    nonmonic_context = Suslin._sl3_realization_input_context(
        nonmonic_target;
        selected_variable = X,
        catalog_metadata = (; fixture_id = "issue-236-nonmonic-local-form"),
    )
    nonmonic_error =
        _captured_error(() -> Suslin._select_sl3_local_form_witness(nonmonic_context))
    @test nonmonic_error isa ArgumentError
    @test occursin("local-form p is not monic", sprint(showerror, nonmonic_error))

    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        multivariate_context;
        local_form_witness = (; entries = (; p = p + one(R), q, r, s)),
    )

    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        variable_context;
        variable_change_metadata = merge(
            variable_change_metadata,
            (; source_matrix = identity_matrix(parent(source_matrix[1, 1]), 3),),
        ),
    )
    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        variable_context;
        variable_change_metadata = merge(
            variable_change_metadata,
            (; selected_variable = X + one(R),),
        ),
    )
    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        normality_shell_context;
        normality_conjugation_metadata = merge(
            normality_shell_context.normality_conjugation_metadata,
            (; source_matrix = identity_matrix(parent(source_matrix[1, 1]), 3),),
        ),
    )

    missing_payload_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-missing-variable-payload"),
        variable_change_metadata = (;
            replay_id = "missing-variable-change-payload",
            source_matrix,
            selected_variable = staged_entry.selected_variable.generator,
            local_form_matrix = multivariate_target,
        ),
    )
    missing_payload_selection =
        Suslin._select_sl3_local_form_witness(missing_payload_context)
    @test missing_payload_selection.support_status == :staged
    @test missing_payload_selection.variable_change_status == :recorded
    @test Suslin._verify_sl3_realization_input_context(missing_payload_context)

    missing_source_matrix_metadata = (;
        replay_id = "missing-source-matrix-replay",
        selected_variable = staged_entry.selected_variable.generator,
        local_form_matrix = multivariate_target,
        replay_steps = ((; kind = :replay_without_source_matrix),),
    )
    missing_source_matrix_status, missing_source_matrix_entries, missing_source_matrix, _ =
        Suslin._sl3_local_witness_metadata_status(
            variable_context,
            missing_source_matrix_metadata,
            staged_entry.selected_variable.generator,
            staged_entry.selected_variable.index,
            "variable-change metadata",
        )
    @test missing_source_matrix_status == :recorded
    @test missing_source_matrix_entries == (; p, q, r, s)
    @test missing_source_matrix == multivariate_target

    missing_source_matrix_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-missing-source-matrix-replay"),
        variable_change_metadata = missing_source_matrix_metadata,
    )
    missing_source_matrix_selection =
        Suslin._select_sl3_local_form_witness(missing_source_matrix_context)
    @test missing_source_matrix_selection.support_status == :staged
    @test missing_source_matrix_selection.variable_change_status == :recorded

    missing_selected_variable_metadata = (;
        replay_id = "missing-selected-variable-replay",
        source_matrix,
        local_form_matrix = multivariate_target,
        replay_steps = ((; kind = :replay_without_selected_variable),),
    )
    missing_selected_variable_status, missing_selected_variable_entries, _, _ =
        Suslin._sl3_local_witness_metadata_status(
            variable_context,
            missing_selected_variable_metadata,
            staged_entry.selected_variable.generator,
            staged_entry.selected_variable.index,
            "variable-change metadata",
        )
    @test missing_selected_variable_status == :recorded
    @test missing_selected_variable_entries == (; p, q, r, s)
    missing_selected_variable_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-missing-selected-variable-replay"),
        variable_change_metadata = missing_selected_variable_metadata,
    )
    missing_selected_variable_selection =
        Suslin._select_sl3_local_form_witness(missing_selected_variable_context)
    @test missing_selected_variable_selection.support_status == :staged
    @test missing_selected_variable_selection.variable_change_status == :recorded

    @test _selection_corruption_is_rejected(
        multivariate_selection,
        (; entries = (; p = p + one(R), q, r, s)),
    )
    @test _selection_corruption_is_rejected(
        multivariate_selection,
        (; monicity_witness = merge(
            multivariate_selection.monicity_witness,
            (; is_monic = false,),
        )),
    )
    @test _selection_corruption_is_rejected(
        variable_selection,
        (; variable_change_status = :missing,),
    )
end
