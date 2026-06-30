module ECPMainlineFixtureCatalog

using Oscar
using Suslin

const ECP_COLUMN_CATALOG_PATH = joinpath(@__DIR__, "ecp_column_cases.jl")
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"

if !isdefined(@__MODULE__, :ECPColumnFixtureCatalog)
    include(ECP_COLUMN_CATALOG_PATH)
end

function _base_cases_by_id()
    return ECPColumnFixtureCatalog.cases_by_id()
end

function _selected_variable(name, generator, index::Int; status = :passes)
    return (;
        name = String(name),
        generator = generator,
        index = index,
        status = status,
    )
end

function _column(entry)
    return tuple((getproperty(entry.column_entries, name) for name in entry.column_order)...)
end

function _column(entry::NamedTuple)
    return tuple((getproperty(entry.entries, name) for name in entry.column_order)...)
end

function _factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _unimodularity(column, R)
    certificate = Suslin.ecp_column_reduction_certificate(collect(column), R)
    product = _factor_product(certificate.factors, R, length(column))
    coefficients = ntuple(idx -> product[length(column), idx], length(column))
    return (;
        status = :replayed,
        witness = certificate,
        coefficients = coefficients,
    )
end

function _gf2_supplied_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    column = _column(entry)
    G = y * column[2] + column[3]
    return (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((;
            probe_id = :gf2_fixture_probe,
            G = G,
            lifted_tail_coefficients = (y, one(R)),
            tilde_G = G,
        ),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
end

function _qq_supplied_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    return (;
        source = :supplied_link_witness,
        residue_probes = (
            (; id = :qq_y_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),
            (; id = :qq_x_probe, kind = :deterministic_fixture, maximal_ideal_generators = (x,)),
        ),
        tail_reductions = (
            (; probe_id = :qq_y_probe, G = y, lifted_tail_coefficients = (zero(R), one(R)), tilde_G = y),
            (; probe_id = :qq_x_probe, G = x, lifted_tail_coefficients = (one(R), zero(R)), tilde_G = x),
        ),
        resultants = (y^2, y + one(R)),
        bezout_coefficients = (
            (; f = zero(R), h = y),
            (; f = one(R), h = -x),
        ),
        coverage_multipliers = (one(R), one(R) - y),
        path_points = (zero(R), y^2 * x, x),
    )
end

function _link_witness_record(entry; variable_order = entry.ring.generators, selected_variable = variable_order[1], supplied_link_witness)
    return Suslin.ecp_link_witness(
        collect(_column(entry)),
        entry.ring.object;
        variable_order = variable_order,
        selected_variable = selected_variable,
        supplied_link_witness = supplied_link_witness,
    )
end

function _link_step_certificate(entry; variable_order = entry.ring.generators, selected_variable = variable_order[1], supplied_link_witness)
    return Suslin.ecp_link_step_certificate(
        collect(_column(entry)),
        entry.ring.object;
        variable_order = variable_order,
        selected_variable = selected_variable,
        supplied_link_witness = supplied_link_witness,
    )
end

function _replace_record_field(record, field::Symbol, value)
    fields = fieldnames(typeof(record))
    values = [getfield(record, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown record field: $(field)")
    values[idx] = value
    return Suslin.ECPLinkWitnessRecord(values...)
end

function _case(;
    id,
    role,
    expected_status,
    ring_constructor,
    ring,
    column_entries,
    column_order,
    selected_variable,
    monicity,
    unimodularity,
    support_evidence,
    source_refs,
    consumer_issue_ids,
    missing_evidence = nothing,
    extras = NamedTuple(),
)
    entry = (;
        id = id,
        role = role,
        expected_status = expected_status,
        ring_constructor = ring_constructor,
        ring = ring,
        column_entries = column_entries,
        column_order = column_order,
        selected_variable = selected_variable,
        monicity = monicity,
        unimodularity = unimodularity,
        support_evidence = support_evidence,
        source_refs = source_refs,
        consumer_issue_ids = consumer_issue_ids,
    )
    missing_evidence === nothing || (entry = merge(entry, (; missing_evidence = missing_evidence)))
    return merge(entry, extras)
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (; id = id, base_case_id = base_case_id, reason = reason))
end

function catalog()
    base_cases = _base_cases_by_id()

    gf2_base = base_cases["ecp-variable-change-monic-gf2"]
    gf2_witness = _gf2_supplied_link_witness(gf2_base)
    gf2_link_record = _link_witness_record(gf2_base; supplied_link_witness = gf2_witness)
    gf2_link_step = _link_step_certificate(gf2_base; supplied_link_witness = gf2_witness)
    gf2_column = _column(gf2_base)
    gf2_unimodularity = _unimodularity(gf2_column, gf2_base.ring.object)
    gf2_lower = Suslin.ecp_column_reduction_certificate(collect(gf2_column), gf2_base.ring.object)

    gf2_hard_slice_case = _case(
        id = "ecp-mainline-gf2-hard-slice",
        role = :gf2_hard_slice,
        expected_status = :supported,
        ring_constructor = gf2_base.ring_constructor,
        ring = gf2_base.ring,
        column_entries = gf2_base.entries,
        column_order = gf2_base.column_order,
        selected_variable = _selected_variable(:x, gf2_base.ring.generators[1], 1),
        monicity = (;
            status = :replayed,
            selected_entry = gf2_base.monicity.selected_entry,
            selected_variable = gf2_base.monicity.variable_name,
            substitution = gf2_base.monicity.substitution,
            transformed_entry = gf2_base.monicity.transformed_entry,
        ),
        unimodularity = gf2_unimodularity,
        support_evidence = (;
            link_witness_status = :replayed,
            link_step_status = :replayed,
            lower_variable_status = :replayed,
            normality_status = :absent,
            sl3_status = :inapplicable,
            link_witness = gf2_link_record,
            link_step = gf2_link_step,
            lower_variable = gf2_lower,
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_4_REF,
            "Provenance: ecp-variable-change-monic-gf2 hard-slice witness family",
        ),
        consumer_issue_ids = ("#185",),
    )

    qq_link_bezout_base = base_cases["ecp-link-bezout-nonunit-witness-qq"]
    qq_link_bezout_column = _column(qq_link_bezout_base)
    qq_link_bezout_case = _case(
        id = "ecp-mainline-qq-link-bezout",
        role = :qq_link_bezout_boundary,
        expected_status = :staged,
        ring_constructor = qq_link_bezout_base.ring_constructor,
        ring = qq_link_bezout_base.ring,
        column_entries = qq_link_bezout_base.entries,
        column_order = qq_link_bezout_base.column_order,
        selected_variable = _selected_variable(:x, qq_link_bezout_base.ring.generators[1], 1),
        monicity = (;
            status = :passes,
            selected_entry = qq_link_bezout_base.monicity.selected_entry,
            selected_variable = qq_link_bezout_base.monicity.variable_name,
            substitution = NamedTuple(),
            transformed_entry = qq_link_bezout_base.monicity.transformed_entry,
        ),
        unimodularity = _unimodularity(qq_link_bezout_column, qq_link_bezout_base.ring.object),
        support_evidence = (;
            link_witness_status = :missing,
            link_step_status = :missing,
            lower_variable_status = :missing,
            normality_status = :absent,
            sl3_status = :absent,
            link_expectation = only(qq_link_bezout_base.witnesses),
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_4_REF,
            "Provenance: ecp-link-bezout-nonunit-witness-qq exact Bezout/resultant witness data",
        ),
        consumer_issue_ids = ("#185", "#88"),
        missing_evidence = (:selected_first_entry_monicity, :link_step, :lower_variable),
    )

    R4, (x4, y4) = Oscar.polynomial_ring(QQ, ["x", "y"])
    length4_entries = (
        a = x4 * y4,
        b = x4 * (one(R4) - y4),
        c = (one(R4) - x4) * y4,
        d = (one(R4) - x4) * (one(R4) - y4),
    )
    length4_column = (length4_entries.a, length4_entries.b, length4_entries.c, length4_entries.d)
    length4_ring = (;
        description = "QQ[x, y]",
        object = R4,
        generator_names = ("x", "y"),
        generators = (x4, y4),
    )
    length4_constructor = (;
        function_name = :polynomial_ring,
        coefficient = "QQ",
        variables = ("x", "y"),
    )
    length4_coupled_case = _case(
        id = "ecp-mainline-length4-coupled-qq",
        role = :length4_all_entry_boundary,
        expected_status = :staged,
        ring_constructor = length4_constructor,
        ring = length4_ring,
        column_entries = length4_entries,
        column_order = (:a, :b, :c, :d),
        selected_variable = _selected_variable(:x, x4, 1),
        monicity = (;
            status = :missing,
            selected_entry = :a,
            selected_variable = :x,
            substitution = NamedTuple(),
            transformed_entry = length4_entries.a,
        ),
        unimodularity = (;
            status = :passes,
            witness = :lagrange_partition_of_unity,
            coefficients = (one(R4), one(R4), one(R4), one(R4)),
        ),
        support_evidence = (;
            link_witness_status = :missing,
            link_step_status = :missing,
            lower_variable_status = :missing,
            normality_status = :absent,
            sl3_status = :absent,
            coupled_support = (;
                all_entries_required = true,
                omitted_corner_points = (
                    (; omitted = :a, common_zero = (zero(R4), zero(R4))),
                    (; omitted = :b, common_zero = (one(R4), zero(R4))),
                    (; omitted = :c, common_zero = (zero(R4), one(R4))),
                    (; omitted = :d, common_zero = (one(R4), one(R4))),
                ),
            ),
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_4_REF,
            "Lagrange-style partition-of-unity boundary requiring all four entries",
        ),
        consumer_issue_ids = ("#185",),
        missing_evidence = (:first_entry_monicity, :link_witness, :link_step, :lower_variable),
    )

    gf2_permuted_base = base_cases["ecp-variable-change-permuted-gf2"]
    gf2_permuted_column = _column(gf2_permuted_base)
    monicity_change_case = _case(
        id = "ecp-mainline-monicity-change-gf2",
        role = :monicity_change_boundary,
        expected_status = :staged,
        ring_constructor = gf2_permuted_base.ring_constructor,
        ring = gf2_permuted_base.ring,
        column_entries = gf2_permuted_base.entries,
        column_order = gf2_permuted_base.column_order,
        selected_variable = _selected_variable(:y, gf2_permuted_base.ring.generators[2], 2),
        monicity = (;
            status = :replayed,
            selected_entry = gf2_permuted_base.monicity.selected_entry,
            selected_variable = gf2_permuted_base.monicity.variable_name,
            substitution = gf2_permuted_base.monicity.substitution,
            transformed_entry = gf2_permuted_base.monicity.transformed_entry,
        ),
        unimodularity = _unimodularity(gf2_permuted_column, gf2_permuted_base.ring.object),
        support_evidence = (;
            link_witness_status = :missing,
            link_step_status = :missing,
            lower_variable_status = :missing,
            normality_status = :absent,
            sl3_status = :inapplicable,
            monicity_provenance = "replayed substitution metadata from ecp-variable-change-permuted-gf2",
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_4_REF,
            "Provenance: ecp-variable-change-permuted-gf2 monicity-changing substitution",
        ),
        consumer_issue_ids = ("#185", "#85"),
        missing_evidence = (:selected_first_entry_monicity, :link_step, :lower_variable),
    )

    qq_sl3_base = base_cases["ecp-monic-first-entry-qq"]
    qq_sl3_witness = _qq_supplied_link_witness(qq_sl3_base)
    qq_sl3_record = _link_witness_record(qq_sl3_base; supplied_link_witness = qq_sl3_witness)
    qq_sl3_column = _column(qq_sl3_base)
    sl3_route_case = _case(
        id = "ecp-mainline-sl3-route-qq",
        role = :sl3_route_boundary,
        expected_status = :staged,
        ring_constructor = qq_sl3_base.ring_constructor,
        ring = qq_sl3_base.ring,
        column_entries = qq_sl3_base.entries,
        column_order = qq_sl3_base.column_order,
        selected_variable = _selected_variable(:x, qq_sl3_base.ring.generators[1], 1),
        monicity = (;
            status = :passes,
            selected_entry = qq_sl3_base.monicity.selected_entry,
            selected_variable = qq_sl3_base.monicity.variable_name,
            substitution = NamedTuple(),
            transformed_entry = qq_sl3_base.monicity.transformed_entry,
        ),
        unimodularity = _unimodularity(qq_sl3_column, qq_sl3_base.ring.object),
        support_evidence = (;
            link_witness_status = :replayed,
            link_step_status = :missing,
            lower_variable_status = :missing,
            normality_status = :missing,
            sl3_status = :missing,
            link_witness = qq_sl3_record,
            sl3_expectation = (;
                route_issue = "#184",
                boundary = :link_realization_before_sl3_driver,
            ),
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_4_REF,
            "Provenance: ecp-monic-first-entry-qq staged toward #184 SL_3 link realization",
        ),
        consumer_issue_ids = ("#185", "#184"),
        missing_evidence = (:link_step, :lower_variable, :normality, :sl3_realization),
    )

    negative_non_unimodular = _negative_control(
        "ecp-mainline-negative-non-unimodular",
        length4_coupled_case.id,
        "break the unimodularity identity while leaving the rest of the metadata intact",
        merge(length4_coupled_case, (;
            unimodularity = merge(length4_coupled_case.unimodularity, (;
                coefficients = (zero(R4), zero(R4), zero(R4), zero(R4)),
            )),
        )),
    )

    bad_path_points = (zero(gf2_base.ring.object), gf2_base.ring.generators[1] + one(gf2_base.ring.object))
    negative_corrupt_link_witness = _negative_control(
        "ecp-mainline-negative-corrupt-link-witness",
        gf2_hard_slice_case.id,
        "mutate the supplied link witness path points so replay fails",
        merge(gf2_hard_slice_case, (;
            support_evidence = merge(gf2_hard_slice_case.support_evidence, (;
                link_witness = _replace_record_field(gf2_link_record, :path_points, bad_path_points),
            )),
        )),
    )

    negative_selected_variable = _negative_control(
        "ecp-mainline-negative-selected-variable-not-generator",
        gf2_hard_slice_case.id,
        "selected variable generator must match the fixture ring generator at its index",
        merge(gf2_hard_slice_case, (;
            selected_variable = _selected_variable(:x, gf2_base.ring.generators[1] + gf2_base.ring.generators[2], 1),
        )),
    )

    negative_supported_without_evidence = _negative_control(
        "ecp-mainline-negative-supported-without-evidence",
        sl3_route_case.id,
        "supported entries must supply replayable link-step and lower-variable evidence",
        merge(sl3_route_case, (;
            expected_status = :supported,
            support_evidence = merge(sl3_route_case.support_evidence, (;
                link_witness_status = :absent,
                link_step_status = :missing,
                lower_variable_status = :missing,
                normality_status = :absent,
                sl3_status = :absent,
            )),
        )),
    )

    return (;
        cases = [
            gf2_hard_slice_case,
            qq_link_bezout_case,
            length4_coupled_case,
            monicity_change_case,
            sl3_route_case,
        ],
        negative_controls = [
            negative_non_unimodular,
            negative_corrupt_link_witness,
            negative_selected_variable,
            negative_supported_without_evidence,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
