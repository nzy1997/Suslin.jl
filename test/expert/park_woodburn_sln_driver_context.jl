using Test
using Oscar
using Suslin

const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")

function _sln_ctx_as_namedtuple(context)
    names = propertynames(context)
    return NamedTuple{names}(Tuple(getproperty(context, name) for name in names))
end

function _sln_ctx_replace_field(context, field::Symbol, value)
    fields = fieldnames(typeof(context))
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown SLn context field $(field)")
    values = [getfield(context, name) for name in fields]
    values[idx] = value
    return Suslin.SLnRecursiveDriverInputContext(values...)
end

function _sln_context_from_entry(entry; variable_order = entry.ring.generators)
    return Suslin._sln_recursive_driver_input_context(
        entry.matrix;
        variable_order = variable_order,
        selected_variable = isempty(entry.ring.generators) ? nothing : entry.ring.generators[1],
        ecp_witness_metadata = entry.peel_steps[1].last_column_ecp,
        final_route_metadata = entry.final_route,
        route_provenance_metadata = entry.route_provenance,
        catalog_id = entry.id,
    )
end

@testset "Park-Woodburn SLn recursive driver context" begin
    if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        include(PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH)
    end
    entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
    catalog = ParkWoodburnSLnDriverFixtureCatalog.catalog()
    negative = Dict(entry.id => entry for entry in catalog.negative_controls)

    mainline = entries["sln-driver-sl4-gf2-ecp-mainline"]
    mainline_ctx = _sln_context_from_entry(mainline)
    @test mainline_ctx.dimension == 4
    @test mainline_ctx.support_classification == :supported
    @test mainline_ctx.staged_reason_code === nothing
    @test mainline_ctx.ecp_evidence_status == :replayed
    @test mainline_ctx.final_route_evidence_status == :replayed
    @test mainline_ctx.route_provenance_status == :recorded
    @test mainline_ctx.last_column == [mainline.matrix[row, 4] for row in 1:4]
    @test mainline_ctx.determinant_status == :one
    @test mainline_ctx.exact_field_status == :supported
    @test Suslin._verify_sln_recursive_driver_input_context(mainline_ctx)

    catalog_only_provenance_ctx = Suslin._sln_recursive_driver_input_context(
        mainline.matrix;
        variable_order = mainline.ring.generators,
        selected_variable = mainline.ring.generators[1],
        ecp_witness_metadata = mainline.peel_steps[1].last_column_ecp,
        final_route_metadata = mainline.final_route,
        catalog_id = mainline.id,
    )
    @test catalog_only_provenance_ctx.route_provenance_status == :missing
    @test Suslin._verify_sln_recursive_driver_input_context(catalog_only_provenance_ctx)

    empty_ecp_shell_ctx = Suslin._sln_recursive_driver_input_context(
        mainline.matrix;
        variable_order = mainline.ring.generators,
        selected_variable = mainline.ring.generators[1],
        ecp_witness_metadata = (;),
        final_route_metadata = mainline.final_route,
        route_provenance_metadata = mainline.route_provenance,
        catalog_id = mainline.id,
    )
    @test empty_ecp_shell_ctx.ecp_evidence_status == :missing
    @test empty_ecp_shell_ctx.staged_reason_code == :missing_ecp_evidence

    recorded_symbol_ecp_ctx = Suslin._sln_recursive_driver_input_context(
        mainline.matrix;
        variable_order = mainline.ring.generators,
        selected_variable = mainline.ring.generators[1],
        ecp_witness_metadata = (; source_case_id = :symbolic_ecp_shell),
        final_route_metadata = mainline.final_route,
        route_provenance_metadata = mainline.route_provenance,
        catalog_id = mainline.id,
    )
    @test recorded_symbol_ecp_ctx.ecp_evidence_status == :recorded
    @test recorded_symbol_ecp_ctx.staged_reason_code == :missing_ecp_evidence

    recorded_tuple_payload_ecp_ctx = Suslin._sln_recursive_driver_input_context(
        mainline.matrix;
        variable_order = mainline.ring.generators,
        selected_variable = mainline.ring.generators[1],
        ecp_witness_metadata = (; replay_steps = ((; kind = :recorded_shell),)),
        final_route_metadata = mainline.final_route,
        route_provenance_metadata = mainline.route_provenance,
        catalog_id = mainline.id,
    )
    @test recorded_tuple_payload_ecp_ctx.ecp_evidence_status == :recorded

    recorded_scalar_payload_ecp_ctx = Suslin._sln_recursive_driver_input_context(
        mainline.matrix;
        variable_order = mainline.ring.generators,
        selected_variable = mainline.ring.generators[1],
        ecp_witness_metadata = (; replay_payload = :recorded_scalar_payload),
        final_route_metadata = mainline.final_route,
        route_provenance_metadata = mainline.route_provenance,
        catalog_id = mainline.id,
    )
    @test recorded_scalar_payload_ecp_ctx.ecp_evidence_status == :recorded

    multistep = entries["sln-driver-sl5-gf2-two-step"]
    multistep_ctx = _sln_context_from_entry(multistep)
    @test multistep_ctx.dimension == 5
    @test multistep_ctx.last_column == [multistep.matrix[row, 5] for row in 1:5]
    @test multistep_ctx.support_classification == :supported
    @test Suslin._verify_sln_recursive_driver_input_context(multistep_ctx)

    legacy = entries["sln-driver-legacy-recursive-column-peel-qq"]
    legacy_ctx = _sln_context_from_entry(legacy)
    @test legacy_ctx.support_classification == :staged
    @test legacy_ctx.staged_reason_code == :missing_ecp_evidence
    @test legacy_ctx.ecp_evidence_status == :missing
    @test Suslin._verify_sln_recursive_driver_input_context(legacy_ctx)

    staged = entries["sln-driver-staged-missing-final-sl3-qq"]
    staged_ctx = _sln_context_from_entry(staged)
    @test staged_ctx.support_classification == :staged
    @test staged_ctx.staged_reason_code == :missing_final_sl3_route
    @test staged_ctx.ecp_evidence_status == :replayed
    @test staged_ctx.final_route_evidence_status == :missing
    @test Suslin._verify_sln_recursive_driver_input_context(staged_ctx)

    unsupported = negative["sln-driver-negative-unsupported-coefficient-ring"]
    unsupported_ctx = Suslin._sln_recursive_driver_input_context(
        identity_matrix(unsupported.ring.object, 4);
        variable_order = unsupported.ring.generators,
        selected_variable = unsupported.ring.generators[1],
        route_provenance_metadata = unsupported.route_provenance,
        catalog_id = unsupported.id,
    )
    @test unsupported_ctx.support_classification == :staged
    @test unsupported_ctx.staged_reason_code == :unsupported_coefficient_ring
    @test unsupported_ctx.exact_field_status == :unsupported
    @test Suslin._verify_sln_recursive_driver_input_context(unsupported_ctx)

    missing_variable_ctx = Suslin._sln_recursive_driver_input_context(
        staged.matrix;
        variable_order = nothing,
        ecp_witness_metadata = staged.peel_steps[1].last_column_ecp,
        final_route_metadata = staged.final_route,
        route_provenance_metadata = staged.route_provenance,
        catalog_id = staged.id,
    )
    @test missing_variable_ctx.staged_reason_code == :missing_variable_metadata

    invalid_order_ctx = Suslin._sln_recursive_driver_input_context(
        staged.matrix;
        variable_order = (staged.ring.generators[1] + one(staged.ring.object),),
        selected_variable = staged.ring.generators[1],
        ecp_witness_metadata = staged.peel_steps[1].last_column_ecp,
        final_route_metadata = staged.final_route,
        route_provenance_metadata = staged.route_provenance,
        catalog_id = staged.id,
    )
    @test invalid_order_ctx.variable_order_status == :missing
    @test invalid_order_ctx.staged_reason_code == :missing_variable_metadata

    invalid_selected_ctx = Suslin._sln_recursive_driver_input_context(
        staged.matrix;
        variable_order = staged.ring.generators,
        selected_variable = staged.ring.generators[1] + one(staged.ring.object),
        ecp_witness_metadata = staged.peel_steps[1].last_column_ecp,
        final_route_metadata = staged.final_route,
        route_provenance_metadata = staged.route_provenance,
        catalog_id = staged.id,
    )
    @test invalid_selected_ctx.selected_variable_status == :missing
    @test invalid_selected_ctx.staged_reason_code == :missing_variable_metadata

    malformed_final_route_ctx = Suslin._sln_recursive_driver_input_context(
        staged.matrix;
        variable_order = staged.ring.generators,
        selected_variable = staged.ring.generators[1],
        ecp_witness_metadata = staged.peel_steps[1].last_column_ecp,
        final_route_metadata = (; status = :replayed, case_id = "malformed-final-route", matrix = (;)),
        route_provenance_metadata = staged.route_provenance,
        catalog_id = staged.id,
    )
    @test malformed_final_route_ctx.final_route_evidence_status == :recorded
    @test malformed_final_route_ctx.staged_reason_code == :missing_final_sl3_route

    det_bad = negative["sln-driver-negative-det-not-one"]
    det_bad_ctx = Suslin._sln_recursive_driver_input_context(
        det_bad.matrix;
        variable_order = det_bad.ring.generators,
        ecp_witness_metadata = det_bad.peel_steps[1].last_column_ecp,
        final_route_metadata = det_bad.final_route,
        route_provenance_metadata = det_bad.route_provenance,
        catalog_id = det_bad.id,
    )
    @test det_bad_ctx.staged_reason_code == :determinant_not_one

    for (field, value) in (
        (:determinant_status, :not_one),
        (:ring_profile, :tampered),
        (:generators, reverse(mainline_ctx.generators)),
        (:last_column, reverse(mainline_ctx.last_column)),
        (:route_provenance_metadata, merge(mainline_ctx.route_provenance_metadata, (; source = "tampered"))),
        (:staged_reason_code, :missing_ecp_evidence),
        (:staged_diagnostic, merge(mainline_ctx.staged_diagnostic, (; message = "tampered"))),
        (:verification, merge(mainline_ctx.verification, (; determinant_status_ok = false))),
    )
        @test !Suslin._verify_sln_recursive_driver_input_context(
            _sln_ctx_replace_field(mainline_ctx, field, value),
        )
    end

    @test_throws ArgumentError Suslin._sln_recursive_driver_input_context(
        identity_matrix(mainline.ring.object, 2),
    )
    @test !Suslin._verify_sln_recursive_driver_input_context((;))
end
