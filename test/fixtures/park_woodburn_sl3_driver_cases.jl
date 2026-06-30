module ParkWoodburnSL3DriverFixtureCatalog

using Oscar
using Suslin

const QUILLEN_MAINLINE_CATALOG_PATH = joinpath(@__DIR__, "quillen_mainline_cases.jl")
const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "quillen_patch_cases.jl")
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"

function _ring_metadata(description, R, generator_names, generators)
    return (;
        description = description,
        object = R,
        generator_names = generator_names,
        generators = generators,
    )
end

function _ordinary_ring_constructor(coefficient, variables)
    return (;
        function_name = :polynomial_ring,
        coefficient = coefficient,
        variables = variables,
    )
end

function _selected_variable(name, generator, index::Int; status = :passes)
    return (;
        name = name,
        generator = generator,
        index = index,
        status = status,
    )
end

function _case(;
    id,
    role,
    expected_status,
    ring_constructor,
    ring,
    matrix,
    selected_variable,
    local_form_status,
    selected_variable_status,
    supplied_witness_status,
    upstream_evidence_status,
    source_refs,
    consumer_issue_ids,
    local_form_witness = nothing,
    supplied_witness = nothing,
    upstream_evidence = nothing,
    staged_reason = nothing,
)
    entry = (;
        id = id,
        role = role,
        expected_status = expected_status,
        ring_constructor = ring_constructor,
        ring = ring,
        matrix = matrix,
        selected_variable = selected_variable,
        local_form_status = local_form_status,
        selected_variable_status = selected_variable_status,
        supplied_witness_status = supplied_witness_status,
        upstream_evidence_status = upstream_evidence_status,
        source_refs = source_refs,
        consumer_issue_ids = consumer_issue_ids,
    )
    local_form_witness === nothing || (entry = merge(entry, (; local_form_witness = local_form_witness)))
    supplied_witness === nothing || (entry = merge(entry, (; supplied_witness = supplied_witness)))
    upstream_evidence === nothing || (entry = merge(entry, (; upstream_evidence = upstream_evidence)))
    staged_reason === nothing || (entry = merge(entry, (; staged_reason = staged_reason)))
    return entry
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (;
        id = id,
        base_case_id = base_case_id,
        reason = reason,
    ))
end

function _mainline_cases_by_id()
    if !isdefined(Main, :QuillenMainlineFixtureCatalog)
        Base.include(Main, QUILLEN_MAINLINE_CATALOG_PATH)
    end
    catalog_module = getfield(Main, :QuillenMainlineFixtureCatalog)
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function _patch_cases_by_id()
    if !isdefined(Main, :QuillenPatchFixtureCatalog)
        Base.include(Main, QUILLEN_PATCH_CATALOG_PATH)
    end
    catalog_module = getfield(Main, :QuillenPatchFixtureCatalog)
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function catalog()
    mainline_cases = _mainline_cases_by_id()
    patch_cases = _patch_cases_by_id()

    RQ, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    RM, (MX, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    RZ, (ZX,) = Oscar.polynomial_ring(ZZ, ["X"])

    qq_univariate_constructor = _ordinary_ring_constructor("QQ", ("X",))
    qq_multivariate_constructor = _ordinary_ring_constructor("QQ", ("X", "r", "g"))
    zz_constructor = _ordinary_ring_constructor("ZZ", ("X",))

    qq_univariate_ring = _ring_metadata("QQ[X]", RQ, ("X",), (X,))
    qq_multivariate_ring = _ring_metadata("QQ[X, r, g]", RM, ("X", "r", "g"), (MX, r, g))
    zz_ring = _ring_metadata("ZZ[X]", RZ, ("X",), (ZX,))

    fast_local_entry = X + one(RQ)
    fast_local_case = _case(
        id = "sl3-driver-univariate-fast-local-qq",
        role = :fast_local_univariate,
        expected_status = :supported,
        ring_constructor = qq_univariate_constructor,
        ring = qq_univariate_ring,
        matrix = elementary_matrix(3, 1, 2, fast_local_entry, RQ),
        selected_variable = _selected_variable("X", X, 1),
        local_form_status = :passes,
        selected_variable_status = :passes,
        supplied_witness_status = :absent,
        upstream_evidence_status = :absent,
        local_form_witness = (;
            entry = fast_local_entry,
            monic_entry_position = 1,
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 184 univariate fast-local driver acceptance shell",
        ),
        consumer_issue_ids = ("#184",),
    )

    multivariate_entry = MX + r * g + one(RM)
    multivariate_case = _case(
        id = "sl3-driver-multivariate-monic-special-form-qq",
        role = :monic_special_form,
        expected_status = :supported,
        ring_constructor = qq_multivariate_constructor,
        ring = qq_multivariate_ring,
        matrix = elementary_matrix(3, 2, 3, multivariate_entry, RM),
        selected_variable = _selected_variable("X", MX, 1),
        local_form_status = :replayed,
        selected_variable_status = :passes,
        supplied_witness_status = :inapplicable,
        upstream_evidence_status = :absent,
        local_form_witness = (;
            polynomial = multivariate_entry,
            monic_entry_position = 2,
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 184 multivariate monic-special-form driver acceptance shell",
        ),
        consumer_issue_ids = ("#184",),
    )

    constructive_mainline = mainline_cases["quillen-constructive-acceptance-gf2"]
    constructive_patch = patch_cases["quillen-constructive-acceptance-gf2"]
    gf2_mainline_case = _case(
        id = "sl3-driver-quillen-mainline-evidence-gf2",
        role = :quillen_mainline_replay,
        expected_status = :supported,
        ring_constructor = constructive_mainline.patch_case.ring_constructor,
        ring = constructive_mainline.patch_case.ring,
        matrix = constructive_mainline.expected_global_product,
        selected_variable = _selected_variable(
            constructive_mainline.patch_case.ring.generator_names[1],
            constructive_mainline.patch_case.ring.generators[1],
            1,
        ),
        local_form_status = :absent,
        selected_variable_status = :passes,
        supplied_witness_status = :absent,
        upstream_evidence_status = :replayed,
        upstream_evidence = (;
            mainline_case_id = constructive_mainline.id,
            patch_case_id = constructive_patch.id,
        ),
        source_refs = (
            PARK_WOODBURN_SECTION_3_REF,
            "Issue 184 GF(2) Quillen mainline replay driver case",
        ),
        consumer_issue_ids = ("#184", "#211"),
    )

    patched_mainline = mainline_cases["quillen-patched-substitution-witness-qq"]
    patched_patch = patch_cases["quillen-patched-substitution-witness-qq"]
    legacy_quillen_case = _case(
        id = "sl3-driver-legacy-quillen-patched-substitution-qq",
        role = :legacy_quillen_patch_replay,
        expected_status = :staged,
        ring_constructor = patched_mainline.patch_case.ring_constructor,
        ring = patched_mainline.patch_case.ring,
        matrix = patched_mainline.expected_global_product,
        selected_variable = _selected_variable(
            patched_mainline.patch_case.ring.generator_names[1],
            patched_mainline.patch_case.ring.generators[1],
            1,
        ),
        local_form_status = :absent,
        selected_variable_status = :passes,
        supplied_witness_status = :absent,
        upstream_evidence_status = :replayed,
        upstream_evidence = (;
            case_id = patched_mainline.id,
            patch_case_id = patched_patch.id,
        ),
        staged_reason = "legacy patched-substitution coverage is recorded for #184 without treating fixture-id equality as the support boundary",
        source_refs = (
            PARK_WOODBURN_SECTION_3_REF,
            "Issue 184 legacy patched-substitution replay driver case",
        ),
        consumer_issue_ids = ("#184", "#214"),
    )

    staged_matrix = elementary_matrix(3, 1, 3, r * g + one(RM), RM)
    staged_case = _case(
        id = "sl3-driver-det-one-no-witness-staged-qq",
        role = :determinant_one_staging_boundary,
        expected_status = :staged,
        ring_constructor = qq_multivariate_constructor,
        ring = qq_multivariate_ring,
        matrix = staged_matrix,
        selected_variable = _selected_variable("X", MX, 1; status = :passes),
        local_form_status = :missing,
        selected_variable_status = :passes,
        supplied_witness_status = :missing,
        upstream_evidence_status = :missing,
        staged_reason = "determinant-one staging case without a recorded local-form or upstream replay witness",
        source_refs = (
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 184 staged determinant-one boundary case",
        ),
        consumer_issue_ids = ("#184",),
    )

    det_not_one_control = _negative_control(
        "sl3-driver-negative-det-not-one",
        fast_local_case.id,
        "matrix determinant is not one",
        merge(fast_local_case, (;
            matrix = matrix(RQ, [
                one(RQ) + X  fast_local_entry  zero(RQ);
                zero(RQ)     one(RQ)          zero(RQ);
                zero(RQ)     zero(RQ)         one(RQ)
            ]),
        )),
    )

    unsupported_ring_control = _negative_control(
        "sl3-driver-negative-unsupported-coefficient-ring",
        fast_local_case.id,
        "coefficient ring is not a field",
        _case(
            id = fast_local_case.id,
            role = :unsupported_coefficient_ring,
            expected_status = :supported,
            ring_constructor = zz_constructor,
            ring = zz_ring,
            matrix = elementary_matrix(3, 1, 2, ZX + one(RZ), RZ),
            selected_variable = _selected_variable("X", ZX, 1),
            local_form_status = :passes,
            selected_variable_status = :passes,
            supplied_witness_status = :absent,
            upstream_evidence_status = :absent,
            local_form_witness = (; p = ZX + one(RZ)),
            source_refs = (
                PARK_WOODBURN_SECTION_5_REF,
                "Issue 184 negative control for unsupported coefficient ring",
            ),
            consumer_issue_ids = ("#184",),
        ),
    )

    bad_selected_variable_control = _negative_control(
        "sl3-driver-negative-selected-variable-not-generator",
        multivariate_case.id,
        "selected variable metadata points at a non-generator polynomial",
        merge(multivariate_case, (;
            selected_variable = _selected_variable("X", MX + one(RM), 1),
        )),
    )

    missing_local_evidence_control = _negative_control(
        "sl3-driver-negative-claimed-local-evidence-missing",
        multivariate_case.id,
        "local-form support is claimed without local_form_witness metadata",
        NamedTuple{Tuple(k for k in keys(multivariate_case) if k != :local_form_witness)}(
            Tuple(v for (k, v) in pairs(multivariate_case) if k != :local_form_witness),
        ),
    )

    supported_without_witness_control = _negative_control(
        "sl3-driver-negative-supported-without-witness",
        staged_case.id,
        "supported status is declared without replayable local or upstream evidence",
        merge(staged_case, (;
            expected_status = :supported,
        )),
    )

    return (;
        cases = [
            fast_local_case,
            multivariate_case,
            gf2_mainline_case,
            legacy_quillen_case,
            staged_case,
        ],
        negative_controls = [
            det_not_one_control,
            unsupported_ring_control,
            bad_selected_variable_control,
            missing_local_evidence_control,
            supported_without_witness_control,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
