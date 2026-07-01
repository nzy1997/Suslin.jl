using Test
using Oscar
using Suslin

const PARK_WOODBURN_SL3_LOCAL_EVIDENCE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

if !isdefined(Main, :SL3MurthyGuptaFixtureCatalog)
    include(PARK_WOODBURN_SL3_LOCAL_EVIDENCE_FIXTURE_PATH)
end

mutable struct Issue237InterruptComparable
    value::Int
end

Base.:(==)(::Issue237InterruptComparable, ::Issue237InterruptComparable) =
    throw(InterruptException())

function _issue237_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}(pair.first => pair.second for pair in kwargs)
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function _issue237_refresh_context_verification(context)
    return _issue237_rebuild(
        context;
        verification = Suslin._sl3_realization_input_context_core_verification(context),
    )
end

function _issue237_refresh_selection_verification(selection)
    return _issue237_rebuild(
        selection;
        verification = Suslin._sl3_local_form_witness_selection_core_verification(selection),
    )
end

function _issue237_context_case()
    R, (X, u, v) = Oscar.polynomial_ring(QQ, ["X", "u", "v"])
    p = X + u * v + one(R)
    q = one(R)
    r = X + u * v
    s = one(R)
    A = matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
    context_metadata = (;
        fixture_id = "issue-237-non-fixture-sl3-context",
        context_issue_id = "#235",
        driver_issue_id = "#184",
        original_matrix_label = :issue237_non_fixture,
    )
    witness_metadata = (;
        entries = (; p, q, r, s),
        source_matrix = A,
        selected_variable = X,
        replay_steps = ((; kind = :issue236_supplied_local_form),),
        witness_issue_id = "#236",
    )
    context = Suslin._sl3_realization_input_context(
        A;
        selected_variable = (; name = "X", generator = X, index = 1, status = :passes),
        catalog_metadata = context_metadata,
        local_form_witness = witness_metadata,
    )
    selection = Suslin._select_sl3_local_form_witness(context)
    return (; R, X, u, v, p, q, r, s, A, context, selection, context_metadata, witness_metadata)
end

function _issue237_variable_change_selection_case()
    case = _issue237_context_case()
    R = case.R
    source_matrix = matrix(R, [
        one(R) zero(R) one(R);
        zero(R) one(R) zero(R);
        zero(R) zero(R) one(R)
    ])
    @assert det(source_matrix) == one(R)
    variable_change_metadata = (;
        replay_id = "issue-237-variable-change-local-form-only",
        source_matrix,
        selected_variable = case.X,
        local_form_matrix = case.A,
        replay_steps = ((; kind = :issue236_supplied_variable_change),),
        witness_issue_id = "#236",
    )
    context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = (;
            name = "X",
            generator = case.X,
            index = 1,
            status = :passes,
        ),
        catalog_metadata = (;
            fixture_id = "issue-237-variable-change-context",
            context_issue_id = "#235",
            driver_issue_id = "#184",
            original_matrix_label = :issue237_variable_change_source,
        ),
        variable_change_metadata,
    )
    selection = Suslin._select_sl3_local_form_witness(context)
    return (; case..., source_matrix, variable_change_metadata, context, selection)
end

function _issue237_nonunit_provider_case()
    catalog = Main.SL3MurthyGuptaFixtureCatalog.catalog()
    entries = Dict(entry.id => entry for entry in catalog.cases)
    entry = entries["mg-q0-nonunit-extracted-bezout-resultant"]
    context_metadata = (;
        fixture_id = entry.id,
        context_issue_id = "#235",
        driver_issue_id = "#184",
        original_matrix_label = :issue237_nonunit_fixture,
    )
    witness_metadata = (;
        entries = entry.entries,
        source_matrix = entry.target,
        selected_variable = entry.variable,
        replay_steps = ((; kind = :issue236_supplied_local_form),),
        witness_issue_id = "#236",
    )
    context = Suslin._sl3_realization_input_context(
        entry.target;
        selected_variable = (;
            name = string(entry.variable),
            generator = entry.variable,
            index = 1,
            status = :passes,
        ),
        catalog_metadata = context_metadata,
        local_form_witness = witness_metadata,
    )
    selection = Suslin._select_sl3_local_form_witness(context)
    provider = Suslin._sl3_murthy_quillen_local_evidence_provider(
        selection;
        metadata = (; provider_test = :fixture_nonunit, route_issue_id = "#184"),
    )
    return (; entry, context, selection, provider, context_metadata, witness_metadata)
end

@testset "Park-Woodburn SL3 Murthy-to-Quillen local evidence provider" begin
    case = _issue237_context_case()
    provider = Suslin._sl3_murthy_quillen_local_evidence_provider(
        case.selection;
        metadata = (; provider_test = :non_fixture, route_issue_id = "#184"),
    )

    @test provider.context == case.context
    @test provider.witness_selection == case.selection
    @test provider.selected_variable == case.X
    @test provider.selected_variable_index == 1
    @test provider.local_product == case.A
    @test provider.murthy_context.target == case.A
    @test Suslin.verify_sl3_local_murthy_input_context(provider.murthy_context)
    @test Suslin.verify_sl3_local_realization(provider.murthy_certificate)
    @test Suslin._verify_murthy_quillen_local_adapter(provider.murthy_adapter)
    @test provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence
    @test provider.staged_diagnostic.status == :supported
    @test provider.denominator_metadata.denominator_product == one(case.R)
    @test provider.denominator_metadata.factor_denominators ==
          [one(case.R) for _ in provider.murthy_adapter.local_factor_replay.factors]
    @test provider.witness_metadata.context_metadata == case.context_metadata
    @test provider.witness_metadata.local_form_witness == case.witness_metadata
    @test provider.witness_metadata.witness_source == :already_special_form
    @test provider.replay_metadata.original_matrix == case.A
    @test provider.replay_metadata.local_product == case.A
    @test provider.replay_metadata.context_metadata == case.context_metadata
    @test length(provider.quillen_local_sequences) == 1
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, provider.quillen_local_sequences)
    @test provider.quillen_local_sequences[1].local_product == case.A
    @test provider.quillen_local_sequences[1].witness_metadata == provider.witness_metadata
    @test Suslin._verify_sl3_murthy_quillen_local_evidence_provider(provider)

    @test_throws ArgumentError Suslin._sl3_realization_input_context(
        case.A;
        selected_variable = (; name = "X", generator = case.X, index = 1, status = :passes),
        catalog_metadata = merge(case.context_metadata, (; expected_status = :unsupported)),
        local_form_witness = case.witness_metadata,
    )

    @test Suslin._same_cached_provenance_data(
        (Ref(1), Ref(2)),
        (Ref(1), Ref(2)),
    )
    @test !Suslin._same_cached_provenance_data(
        (Ref(1),),
        (Ref(1), Ref(2)),
    )
    @test Suslin._same_cached_provenance_data(
        [Ref(1), Ref(2)],
        [Ref(1), Ref(2)],
    )
    @test !Suslin._same_cached_provenance_data(
        [Ref(1)],
        [Ref(1), Ref(2)],
    )
    @test_throws InterruptException Suslin._same_cached_provenance_data(
        Issue237InterruptComparable(1),
        Issue237InterruptComparable(1),
    )

    retagged_witness_metadata = merge(
        case.witness_metadata,
        (;
            witness_issue_id = "#236-forged",
            replay_steps = ((; kind = :forged_issue236_local_form_replay),),
        ),
    )
    retagged_witness_context = _issue237_refresh_context_verification(
        _issue237_rebuild(case.context; local_form_witness = retagged_witness_metadata),
    )
    @test !Suslin._verify_sl3_realization_input_context(retagged_witness_context)

    retagged_witness_selection = _issue237_refresh_selection_verification(
        _issue237_rebuild(case.selection; local_form_witness = retagged_witness_metadata),
    )
    @test_throws ArgumentError Suslin._sl3_murthy_quillen_local_evidence_provider(
        retagged_witness_selection,
    )
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; witness_selection = retagged_witness_selection),
    )

    variable_change_case = _issue237_variable_change_selection_case()
    @test variable_change_case.selection.witness_source == :variable_change
    @test variable_change_case.selection.local_form_matrix != variable_change_case.context.matrix
    @test Suslin._verify_sl3_local_form_witness_selection(variable_change_case.selection)
    @test_throws ArgumentError Suslin._sl3_murthy_quillen_local_evidence_provider(
        variable_change_case.selection,
    )

    no_witness_context = Suslin._sl3_realization_input_context(
        case.A;
        selected_variable = (; name = "X", generator = case.X, index = 1, status = :passes),
        catalog_metadata = merge(
            case.context_metadata,
            (; fixture_id = "issue-237-no-issue236-witness"),
        ),
    )
    no_witness_selection = Suslin._select_sl3_local_form_witness(no_witness_context)
    @test no_witness_selection.support_status == :supported
    @test no_witness_selection.witness_source == :already_special_form
    @test no_witness_selection.local_form_witness === nothing
    @test_throws ArgumentError Suslin._sl3_murthy_quillen_local_evidence_provider(
        no_witness_selection,
    )

    shell_witness_metadata = (;
        source_matrix = case.A,
        selected_variable = case.X,
        replay_steps = ((; kind = :issue236_shell_without_local_form_entries),),
        witness_issue_id = "#236",
    )
    shell_witness_context = Suslin._sl3_realization_input_context(
        case.A;
        selected_variable = (; name = "X", generator = case.X, index = 1, status = :passes),
        catalog_metadata = merge(
            case.context_metadata,
            (; fixture_id = "issue-237-shell-issue236-witness"),
        ),
        local_form_witness = shell_witness_metadata,
    )
    shell_witness_selection = Suslin._select_sl3_local_form_witness(shell_witness_context)
    @test shell_witness_selection.support_status == :supported
    @test shell_witness_selection.witness_source == :already_special_form
    @test shell_witness_selection.local_form_witness == shell_witness_metadata
    @test_throws ArgumentError Suslin._sl3_murthy_quillen_local_evidence_provider(
        shell_witness_selection,
    )

    bad_selection = _issue237_rebuild(case.selection; selected_variable_index = 2)
    @test !Suslin._verify_sl3_local_form_witness_selection(bad_selection)
    @test_throws ArgumentError Suslin._sl3_murthy_quillen_local_evidence_provider(bad_selection)
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; witness_selection = bad_selection),
    )

    bad_certificate = Suslin.SL3LocalRealizationCertificate(
        provider.murthy_certificate.target,
        provider.murthy_certificate.branch,
        reverse(provider.murthy_certificate.factors),
        provider.murthy_certificate.selected_variable,
        provider.murthy_certificate.witness,
    )
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; murthy_certificate = bad_certificate),
    )

    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; selected_variable = case.u),
    )

    bad_denominator_metadata = merge(
        provider.denominator_metadata,
        (; denominator_product = provider.denominator_metadata.denominator_product + case.X),
    )
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; denominator_metadata = bad_denominator_metadata),
    )

    bad_context = _issue237_rebuild(
        case.context;
        catalog_metadata = (; fixture_id = "tampered-context"),
    )
    @test !Suslin._verify_sl3_realization_input_context(bad_context)
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; context = bad_context),
    )

    retagged_metadata = merge(
        case.context_metadata,
        (; fixture_id = "retagged-context", review_issue_id = "#237"),
    )
    retagged_context = _issue237_refresh_context_verification(
        _issue237_rebuild(case.context; catalog_metadata = retagged_metadata),
    )
    @test !Suslin._verify_sl3_realization_input_context(retagged_context)
    forged_identity = (;
        source = :sl3_realization_input_context,
        snapshot = (; malformed = true),
    )
    synchronized_retagged_context = _issue237_refresh_context_verification(
        _issue237_rebuild(
            case.context;
            catalog_metadata = retagged_metadata,
            creation_identity = forged_identity,
        ),
    )
    @test !Suslin._verify_sl3_realization_input_context(synchronized_retagged_context)
    reconstructed_context = Suslin._sl3_realization_input_context(
        case.A;
        selected_variable = (;
            name = "X",
            generator = case.X,
            index = 1,
            status = :passes,
        ),
        catalog_metadata = retagged_metadata,
        local_form_witness = case.witness_metadata,
    )
    @test Suslin._verify_sl3_realization_input_context(reconstructed_context)

    mutable_matrix_context = Suslin._sl3_realization_input_context(
        case.A;
        selected_variable = (;
            name = "X",
            generator = case.X,
            index = 1,
            status = :passes,
        ),
        catalog_metadata = case.context_metadata,
        local_form_witness = case.witness_metadata,
    )
    mutable_matrix_context.matrix[1, 1] += one(case.R)
    mutable_matrix_context.matrix[2, 1] += one(case.R)
    @test det(mutable_matrix_context.matrix) == one(case.R)
    refreshed_matrix_context = _issue237_refresh_context_verification(mutable_matrix_context)
    @test !Suslin._verify_sl3_realization_input_context(refreshed_matrix_context)

    retagged_selection = _issue237_rebuild(case.selection; context = retagged_context)
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(
            provider;
            context = retagged_context,
            witness_selection = retagged_selection,
        ),
    )

    nonunit_case = _issue237_nonunit_provider_case()
    nonunit_provider = nonunit_case.provider
    @test nonunit_provider.murthy_certificate.branch == :murthy_q0_nonunit_bezout_resultant
    original_reduction = nonunit_provider.murthy_certificate.witness.reduction
    tampered_witness_source =
        original_reduction.witness_source == :extracted_bezout_witness ?
        :supplied_bezout_witness :
        :extracted_bezout_witness
    tampered_reduction = _issue237_rebuild(
        original_reduction;
        witness_source = tampered_witness_source,
    )
    tampered_certificate = Suslin.SL3LocalRealizationCertificate(
        nonunit_provider.murthy_certificate.target,
        nonunit_provider.murthy_certificate.branch,
        nonunit_provider.murthy_certificate.factors,
        nonunit_provider.murthy_certificate.selected_variable,
        merge(
            nonunit_provider.murthy_certificate.witness,
            (; reduction = tampered_reduction),
        ),
    )
    @test Suslin.verify_sl3_local_realization(tampered_certificate)
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(nonunit_provider; murthy_certificate = tampered_certificate),
    )
end
