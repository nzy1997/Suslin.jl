module QuillenMainlineFixtureCatalog

using Oscar
using Suslin

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "quillen_patch_cases.jl")
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"

if !isdefined(Main, :QuillenPatchFixtureCatalog)
    Base.include(Main, QUILLEN_PATCH_CATALOG_PATH)
end

function _patch_catalog_module()
    return getfield(Main, :QuillenPatchFixtureCatalog)
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (;
        id = id,
        base_case_id = base_case_id,
        reason = reason,
    ))
end

function _raw_denominator_record(index::Int, data; source_ref)
    return (;
        local_index = index,
        denominator = data.denominator,
        coverage_multiplier = data.coverage_multiplier,
        source_ref = source_ref,
    )
end

function _raw_denominator_provenance(denominator_data; source_ref)
    return (;
        records = Tuple(
            _raw_denominator_record(index, data; source_ref = source_ref)
            for (index, data) in enumerate(denominator_data)
        ),
        source_ref = source_ref,
    )
end

function _denominator_cover_term(; denominator, multiplier)
    return (;
        denominator = denominator,
        multiplier = multiplier,
        term = denominator * multiplier,
    )
end

function _denominator_cover_from_data(denominator_data, R)
    denominators = Tuple(data.denominator for data in denominator_data)
    multipliers = Tuple(data.coverage_multiplier for data in denominator_data)
    coverage_terms = Tuple(_denominator_cover_term(;
        denominator = data.denominator,
        multiplier = data.coverage_multiplier,
    ) for data in denominator_data)
    coverage_sum = foldl(
        (acc, data) -> acc + data.denominator * data.coverage_multiplier,
        denominator_data;
        init = zero(R),
    )
    return (;
        denominators = denominators,
        multipliers = multipliers,
        coverage_terms = coverage_terms,
        coverage_sum = coverage_sum,
    )
end

function _local_evidence_record(; kind, source, sequence_status, factor)
    return (;
        kind = kind,
        source = source,
        sequence_status = sequence_status,
        factor = factor,
    )
end

function _local_evidence(factors, expected_product, source_ref; sequence_status = :recorded)
    records = Tuple(_local_evidence_record(
        kind = :local_factor,
        source = source_ref,
        sequence_status = sequence_status,
        factor = factor,
    ) for factor in factors)
    return (;
        factors = factors,
        expected_product = expected_product,
        records = records,
    )
end

function _patched_substitution_chain_step(input, variable, denominator, exponent_l::Int, multiplier)
    output = input - variable * denominator^exponent_l * multiplier
    return (;
        input = input,
        denominator = denominator,
        exponent_l = exponent_l,
        multiplier = multiplier,
        output = output,
    )
end

function _patched_substitution_chain(
    variable,
    denominator,
    exponent::Int,
    shift,
    expected_matrix,
    coverage_terms,
    exponent_l::Int;
    start = variable,
    status = :complete,
    staging_reason = nothing,
)
    steps = ()
    current = start
    for term in coverage_terms
        step = _patched_substitution_chain_step(
            current,
            variable,
            term.denominator,
            exponent_l,
            term.multiplier,
        )
        steps = (steps..., step)
        current = step.output
    end
    return (;
        variable = variable,
        denominator = denominator,
        exponent = exponent,
        shift = shift,
        expected_matrix = expected_matrix,
        status = status,
        staging_reason = staging_reason,
        sign_convention = :park_woodburn_minus,
        start = start,
        steps = steps,
        final_variable = current,
    )
end

function _base_term_evidence(; status::Symbol, source_ref)
    return (;
        status = status,
        source_ref = source_ref,
    )
end

function _mainline_case(;
    id,
    patch_case,
    denominator_cover,
    raw_denominator_provenance,
    local_evidence,
    patched_substitution_chain,
    base_term_evidence,
    source_refs,
    consumer_issue_ids,
    murthy_adapter_handoff = nothing,
)
    return (;
        id = id,
        source_patch_catalog = :quillen_patch_cases,
        source_patch_case_id = id,
        patch_case = patch_case,
        expected_global_product = patch_case.target_matrix,
        denominator_cover = denominator_cover,
        raw_denominator_provenance = raw_denominator_provenance,
        local_evidence = local_evidence,
        patched_substitution_chain = patched_substitution_chain,
        base_term_evidence = base_term_evidence,
        source_refs = source_refs,
        consumer_issue_ids = consumer_issue_ids,
        murthy_adapter_handoff = murthy_adapter_handoff,
    )
end

function _append_common_issue_refs(patch_case_issue_ids, issue_ids...)
    all_ids = collect(patch_case_issue_ids)
    for id in issue_ids
        id in all_ids || push!(all_ids, id)
    end
    return Tuple(all_ids)
end

function catalog()
    patch_module = _patch_catalog_module()
    patch_entries = Base.invokelatest(getfield(patch_module, :cases_by_id))

    two_open = patch_entries["quillen-two-open-cover-qq"]
    two_open_ring = two_open.ring.object
    two_open_cover = _denominator_cover_from_data(two_open.denominator_data, two_open_ring)
    two_open_local_factors = Tuple(factor.factor for factor in two_open.local_factors)
    two_open_local_evidence = _local_evidence(
        two_open_local_factors,
        two_open.expected.global_correction,
        :two_open_record,
        sequence_status = :recorded,
    )
    two_open_chain = _patched_substitution_chain(
        two_open.substitution_variable,
        two_open.denominator_data[1].denominator,
        1,
        zero(two_open_ring),
        two_open.target_matrix,
        two_open_cover.coverage_terms,
        1,
    )
    two_open_mainline = _mainline_case(
        id = two_open.id,
        patch_case = two_open,
        denominator_cover = two_open_cover,
        raw_denominator_provenance = _raw_denominator_provenance(
            two_open.denominator_data;
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        local_evidence = two_open_local_evidence,
        patched_substitution_chain = two_open_chain,
        base_term_evidence = _base_term_evidence(
            status = :staged,
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        source_refs = (PARK_WOODBURN_SECTION_3_REF, "Issue 99 two-open cover"),
        consumer_issue_ids = _append_common_issue_refs(
            two_open.consumer_issue_ids,
            "#183",
        ),
    )

    nontrivial = patch_entries["quillen-nontrivial-multipliers-qq"]
    nontrivial_ring = nontrivial.ring.object
    nontrivial_cover = _denominator_cover_from_data(nontrivial.denominator_data, nontrivial_ring)
    nontrivial_local_factors = Tuple(factor.factor for factor in nontrivial.local_factors)
    nontrivial_local_evidence = _local_evidence(
        nontrivial_local_factors,
        nontrivial.expected.global_correction,
        :nontrivial_record,
    )
    nontrivial_chain = _patched_substitution_chain(
        nontrivial.substitution_variable,
        nontrivial.denominator_data[1].denominator,
        1,
        zero(nontrivial_ring),
        nontrivial.target_matrix,
        nontrivial_cover.coverage_terms,
        1,
    )
    nontrivial_mainline = _mainline_case(
        id = nontrivial.id,
        patch_case = nontrivial,
        denominator_cover = nontrivial_cover,
        raw_denominator_provenance = _raw_denominator_provenance(
            nontrivial.denominator_data;
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        local_evidence = nontrivial_local_evidence,
        patched_substitution_chain = nontrivial_chain,
        base_term_evidence = _base_term_evidence(
            status = :assumes_identity,
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        source_refs = (PARK_WOODBURN_SECTION_3_REF, "Issue 99 nontrivial coverage multiplier"),
        consumer_issue_ids = _append_common_issue_refs(
            nontrivial.consumer_issue_ids,
            "#183",
        ),
    )

    patched = patch_entries["quillen-patched-substitution-witness-qq"]
    patched_ring = patched.ring.object
    patched_cover = _denominator_cover_from_data(patched.denominator_data, patched_ring)
    patched_local_factors = Tuple(factor.factor for factor in patched.local_factors)
    patched_local_evidence = _local_evidence(
        patched_local_factors,
        patched.expected.global_correction,
        :patched_record,
        sequence_status = :placeholder_for_issue_214,
    )
    patched_chain = _patched_substitution_chain(
        patched.substitution_variable,
        patched.denominator_data[1].denominator,
        patched.patched_substitution_witness.exponent,
        patched.patched_substitution_witness.shift,
        patched.patched_substitution_witness.expected_matrix,
        patched_cover.coverage_terms,
        patched.patched_substitution_witness.exponent,
        status = :staged,
        staging_reason = "powered cover multipliers are staged until #183 supplies exponent-l cover evidence",
    )
    patched_mainline = _mainline_case(
        id = patched.id,
        patch_case = patched,
        denominator_cover = patched_cover,
        raw_denominator_provenance = _raw_denominator_provenance(
            patched.denominator_data;
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        local_evidence = patched_local_evidence,
        patched_substitution_chain = patched_chain,
        base_term_evidence = _base_term_evidence(
            status = :staged,
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_3_REF,
            "Issue 99 patched substitution witness",
            "Issue 214 placeholder sequence metadata",
        ),
        consumer_issue_ids = _append_common_issue_refs(
            patched.consumer_issue_ids,
            "#183",
            "#214",
        ),
    )

    constructive = patch_entries["quillen-constructive-acceptance-gf2"]
    constructive_ring = constructive.ring.object
    constructive_cover = _denominator_cover_from_data(constructive.denominator_data, constructive_ring)
    constructive_local_factors = Tuple(factor.factor for factor in constructive.local_factors)
    constructive_local_evidence = _local_evidence(
        constructive_local_factors,
        constructive.expected.global_correction,
        :constructive_record,
    )
    constructive_chain = _patched_substitution_chain(
        constructive.substitution_variable,
        constructive.denominator_data[1].denominator,
        1,
        zero(constructive_ring),
        constructive.target_matrix,
        constructive_cover.coverage_terms,
        1,
    )
    constructive_mainline = _mainline_case(
        id = constructive.id,
        patch_case = constructive,
        denominator_cover = constructive_cover,
        raw_denominator_provenance = _raw_denominator_provenance(
            constructive.denominator_data;
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        local_evidence = constructive_local_evidence,
        patched_substitution_chain = constructive_chain,
        base_term_evidence = _base_term_evidence(
            status = :assumes_identity,
            source_ref = PARK_WOODBURN_SECTION_3_REF,
        ),
        source_refs = (PARK_WOODBURN_SECTION_3_REF, "Issue 99 constructive acceptance"),
        consumer_issue_ids = _append_common_issue_refs(
            constructive.consumer_issue_ids,
            "#183",
            "#211",
        ),
        murthy_adapter_handoff = (
            status = :staged_until_adapter,
            issue_id = "#211",
            accepts_denominator_one_factors = true,
        ),
    )

    mutated_two_open_factor = two_open_local_factors[1] * elementary_matrix(
        two_open.size,
        1,
        2,
        one(two_open_ring),
        two_open_ring,
    )
    tampered_local_evidence = (
        factors = (
            mutated_two_open_factor,
            two_open_local_factors[2],
        ),
        expected_product = two_open.expected.global_correction,
        records = (
            _local_evidence_record(
                kind = :local_factor,
                source = :two_open_record,
                sequence_status = :recorded,
                factor = mutated_two_open_factor,
            ),
            two_open_mainline.local_evidence.records[2],
        ),
    )

    bad_cover = _negative_control(
        "quillen-mainline-uncovered-denominator-control",
        two_open_mainline.id,
        "mutated coverage multiplier makes denominator cover sum differ from one",
        merge(two_open_mainline, (;
            denominator_cover = _denominator_cover_from_data((
                two_open.denominator_data[1],
                merge(two_open.denominator_data[2], (;
                    coverage_multiplier = two_open.ring.generators[2],
                )),
            ), two_open.ring.object),
        )),
    )
    bad_local = _negative_control(
        "quillen-mainline-tampered-local-evidence-control",
        two_open_mainline.id,
        "mutated first local evidence factor while keeping expected product",
        merge(two_open_mainline, (;
            local_evidence = tampered_local_evidence,
        )),
    )

    return (;
        cases = [
            two_open_mainline,
            nontrivial_mainline,
            patched_mainline,
            constructive_mainline,
        ],
        negative_controls = [bad_cover, bad_local],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
