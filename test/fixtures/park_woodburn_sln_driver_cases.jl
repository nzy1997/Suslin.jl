module ParkWoodburnSLnDriverFixtureCatalog

using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_sl3_driver_cases.jl")
const PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_polynomial_cases.jl")
const PARK_WOODBURN_SLN_REDUCTION_REF =
    "refs/arXiv-alg-geom9405003v1 Section \"Reduction to SL_3(k[x_1,...,x_m])\""
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"

if !isdefined(@__MODULE__, :ParkWoodburnSL3DriverFixtureCatalog)
    include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
end
if !isdefined(@__MODULE__, :ParkWoodburnPolynomialFixtureCatalog)
    include(PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH)
end

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

function _route_provenance(source; route = :park_woodburn_sln_recursive_driver)
    return (;
        route = route,
        reduction_ref = PARK_WOODBURN_SLN_REDUCTION_REF,
        source = source,
    )
end

function _case(;
    id,
    support_role,
    expected_status,
    route_provenance,
    ring_constructor,
    ring,
    matrix,
    expected_peel_count,
    descent_dimensions,
    peel_steps,
    final_route,
    source_refs,
    consumer_issue_ids,
    staged_reason_codes = (),
)
    return (;
        id = id,
        support_role = support_role,
        expected_status = expected_status,
        staged_reason_codes = staged_reason_codes,
        route_provenance = route_provenance,
        ring_constructor = ring_constructor,
        ring = ring,
        matrix = matrix,
        expected_peel_count = expected_peel_count,
        descent_dimensions = descent_dimensions,
        peel_steps = peel_steps,
        final_route = final_route,
        source_refs = source_refs,
        consumer_issue_ids = consumer_issue_ids,
    )
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (;
        id = id,
        base_case_id = base_case_id,
        reason = reason,
    ))
end

function _factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _inverse_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in reverse(collect(factors))
        product *= inv(factor)
    end
    return product
end

function _matrix_column(matrix_value, col::Int)
    return [matrix_value[row, col] for row in 1:nrows(matrix_value)]
end

function _target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _build_peel_step(
    next_block,
    last_column;
    bottom_row_entries,
    ecp_status = :replayed,
    ecp_source_case_id = nothing,
)
    R = base_ring(next_block)
    d = nrows(next_block) + 1
    certificate = ecp_status == :replayed ?
        Suslin.ecp_column_reduction_certificate(collect(last_column), R) :
        nothing
    left_factors = certificate === nothing ? typeof(identity_matrix(R, d))[] : certificate.factors
    after_left = block_embedding(next_block, d, collect(1:(d - 1)))
    for col in 1:(d - 1)
        after_left[d, col] = bottom_row_entries[col]
    end
    input_matrix = _inverse_factor_product(left_factors, R, d) * after_left
    right_factors = typeof(identity_matrix(R, d))[]
    for col in 1:(d - 1)
        coeff = -after_left[d, col]
        coeff == zero(R) || push!(right_factors, elementary_matrix(d, d, col, coeff, R))
    end
    right_product = _factor_product(right_factors, R, d)
    peeled_matrix = after_left * right_product
    return (;
        dimension = d,
        input_matrix = input_matrix,
        last_column = collect(last_column),
        last_column_ecp = (;
            status = ecp_status,
            certificate = certificate,
            source_case_id = ecp_source_case_id,
            target_column = _target_column(R, d),
        ),
        right_clearing = (;
            status = :replayed,
            after_left_matrix = after_left,
            right_factors = right_factors,
            peeled_matrix = peeled_matrix,
            next_block = next_block,
        ),
        next_block = next_block,
    )
end

function _build_legacy_peel_step(input_matrix, next_block)
    R = base_ring(input_matrix)
    d = nrows(input_matrix)
    after_left = block_embedding(next_block, d, collect(1:(d - 1)))
    right_factors = typeof(identity_matrix(R, d))[]
    return (;
        dimension = d,
        input_matrix = input_matrix,
        last_column = _matrix_column(input_matrix, d),
        last_column_ecp = (;
            status = :absent,
            certificate = nothing,
            source_case_id = nothing,
            target_column = _target_column(R, d),
        ),
        right_clearing = (;
            status = :replayed,
            after_left_matrix = after_left,
            right_factors = right_factors,
            peeled_matrix = after_left,
            next_block = next_block,
        ),
        next_block = next_block,
    )
end

function _sl3_cases_by_id()
    return Base.invokelatest(getfield(ParkWoodburnSL3DriverFixtureCatalog, :cases_by_id))
end

function _polynomial_cases_by_id()
    return Base.invokelatest(getfield(ParkWoodburnPolynomialFixtureCatalog, :cases_by_id))
end

function _final_route_replayed(case_id, final_block)
    return (;
        status = :replayed,
        case_id = case_id,
        source_catalog = :park_woodburn_sl3_driver_cases,
        matrix = final_block,
    )
end

function catalog()
    sl3_cases = _sl3_cases_by_id()
    polynomial_cases = _polynomial_cases_by_id()

    legacy_source = polynomial_cases["pw-poly-recursive-column-peel-sl3-qq"]
    legacy_final = polynomial_cases["pw-poly-univariate-sl3-fast-local-qq"]
    legacy_step = _build_legacy_peel_step(legacy_source.matrix, legacy_final.matrix)
    legacy_case = _case(
        id = "sln-driver-legacy-recursive-column-peel-qq",
        support_role = :legacy_regression,
        expected_status = :staged,
        staged_reason_codes = (:legacy_regression_only,),
        route_provenance = _route_provenance(
            "legacy polynomial recursive column-peel fixture $(legacy_source.id)";
            route = :legacy_polynomial_recursive_column_peel,
        ),
        ring_constructor = legacy_source.ring_constructor,
        ring = legacy_source.ring,
        matrix = legacy_source.matrix,
        expected_peel_count = 1,
        descent_dimensions = (4, 3),
        peel_steps = (legacy_step,),
        final_route = (;
            status = :legacy_regression,
            source_catalog = :park_woodburn_polynomial_cases,
            source_case_id = legacy_source.id,
            final_case_id = legacy_final.id,
            matrix = legacy_final.matrix,
        ),
        source_refs = (
            PARK_WOODBURN_SLN_REDUCTION_REF,
            PARK_WOODBURN_SECTION_5_REF,
            legacy_source.source_refs...,
        ),
        consumer_issue_ids = ("#186",),
    )

    gf2_final_case_id = "sl3-driver-quillen-mainline-evidence-gf2"
    gf2_final = sl3_cases[gf2_final_case_id]
    R2 = gf2_final.ring.object
    X2, r2, s2 = gf2_final.ring.generators
    gf2_col4 = (
        X2 + r2^2,
        X2 * r2 + X2 + one(R2),
        X2^2 + X2 * r2 + r2 + one(R2),
        zero(R2),
    )
    gf2_step4 = _build_peel_step(
        gf2_final.matrix,
        gf2_col4;
        bottom_row_entries = (s2, r2, X2 + s2),
        ecp_source_case_id = gf2_final_case_id,
    )
    gf2_sl4_case = _case(
        id = "sln-driver-sl4-gf2-ecp-mainline",
        support_role = :issue186_mainline,
        expected_status = :supported,
        route_provenance = _route_provenance(
            "issue 186 SL4 peel to supported issue 184 GF(2) SL3 route",
        ),
        ring_constructor = gf2_final.ring_constructor,
        ring = gf2_final.ring,
        matrix = gf2_step4.input_matrix,
        expected_peel_count = 1,
        descent_dimensions = (4, 3),
        peel_steps = (gf2_step4,),
        final_route = _final_route_replayed(gf2_final_case_id, gf2_final.matrix),
        source_refs = (
            PARK_WOODBURN_SLN_REDUCTION_REF,
            PARK_WOODBURN_SECTION_4_REF,
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 186 SL4 GF(2) replayed ECP peel to issue 184 SL3 support",
        ),
        consumer_issue_ids = ("#186",),
    )

    gf2_col5 = (
        X2 + r2^2,
        X2 * r2 + X2 + one(R2),
        X2^2 + X2 * r2 + r2 + one(R2),
        zero(R2),
        zero(R2),
    )
    gf2_step5 = _build_peel_step(
        gf2_sl4_case.matrix,
        gf2_col5;
        bottom_row_entries = (one(R2), s2, r2 + s2, X2 * s2 + one(R2)),
        ecp_source_case_id = gf2_sl4_case.id,
    )
    gf2_sl5_case = _case(
        id = "sln-driver-sl5-gf2-two-step",
        support_role = :issue186_mainline,
        expected_status = :supported,
        route_provenance = _route_provenance(
            "issue 186 SL5 two-step peel through the supported GF(2) SL4 driver fixture",
        ),
        ring_constructor = gf2_final.ring_constructor,
        ring = gf2_final.ring,
        matrix = gf2_step5.input_matrix,
        expected_peel_count = 2,
        descent_dimensions = (5, 4, 3),
        peel_steps = (gf2_step5, gf2_step4),
        final_route = _final_route_replayed(gf2_final_case_id, gf2_final.matrix),
        source_refs = (
            PARK_WOODBURN_SLN_REDUCTION_REF,
            PARK_WOODBURN_SECTION_4_REF,
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 186 SL5 GF(2) two-step replayed ECP peel to issue 184 SL3 support",
        ),
        consumer_issue_ids = ("#186",),
    )

    qq_final_case_id = "sl3-driver-multivariate-monic-special-form-qq"
    qq_final = sl3_cases[qq_final_case_id]
    RQ = qq_final.ring.object
    XQ, rQ, gQ = qq_final.ring.generators
    qq_col4 = (
        XQ + rQ,
        XQ * gQ + one(RQ),
        rQ + gQ,
        one(RQ),
    )
    qq_step4 = _build_peel_step(
        qq_final.matrix,
        qq_col4;
        bottom_row_entries = (rQ, gQ, XQ + one(RQ)),
        ecp_source_case_id = qq_final_case_id,
    )
    qq_sl4_case = _case(
        id = "sln-driver-sl4-final-sl3-evidence-qq",
        support_role = :issue186_mainline,
        expected_status = :supported,
        route_provenance = _route_provenance(
            "issue 186 SL4 peel to supported QQ multivariate SL3 driver evidence",
        ),
        ring_constructor = qq_final.ring_constructor,
        ring = qq_final.ring,
        matrix = qq_step4.input_matrix,
        expected_peel_count = 1,
        descent_dimensions = (4, 3),
        peel_steps = (qq_step4,),
        final_route = _final_route_replayed(qq_final_case_id, qq_final.matrix),
        source_refs = (
            PARK_WOODBURN_SLN_REDUCTION_REF,
            PARK_WOODBURN_SECTION_4_REF,
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 186 SL4 QQ replayed ECP peel to issue 184 SL3 monic special-form support",
        ),
        consumer_issue_ids = ("#186",),
    )

    staged_final_block = elementary_matrix(3, 1, 3, XQ * rQ + gQ + one(RQ), RQ)
    staged_col4 = (
        XQ + gQ,
        rQ * gQ + one(RQ),
        XQ + rQ + one(RQ),
        one(RQ),
    )
    staged_step4 = _build_peel_step(
        staged_final_block,
        staged_col4;
        bottom_row_entries = (gQ, XQ, rQ + one(RQ)),
        ecp_source_case_id = "sln-driver-staged-missing-final-sl3-qq",
    )
    staged_missing_case = _case(
        id = "sln-driver-staged-missing-final-sl3-qq",
        support_role = :staged_issue186_candidate,
        expected_status = :staged,
        staged_reason_codes = (:missing_final_sl3_route,),
        route_provenance = _route_provenance(
            "issue 186 determinant-one SL4 staged candidate with replayed ECP but no final SL3 route",
        ),
        ring_constructor = qq_final.ring_constructor,
        ring = qq_final.ring,
        matrix = staged_step4.input_matrix,
        expected_peel_count = 1,
        descent_dimensions = (4, 3),
        peel_steps = (staged_step4,),
        final_route = (;
            status = :missing,
            source_catalog = :park_woodburn_sl3_driver_cases,
            matrix = staged_final_block,
        ),
        source_refs = (
            PARK_WOODBURN_SLN_REDUCTION_REF,
            PARK_WOODBURN_SECTION_4_REF,
            "Issue 186 staged SL4 QQ replayed ECP peel without final SL3 route evidence",
        ),
        consumer_issue_ids = ("#186",),
    )

    RZ, (XZ,) = Oscar.polynomial_ring(ZZ, ["X"])
    zz_ring_constructor = _ordinary_ring_constructor("ZZ", ("X",))
    zz_ring = _ring_metadata("ZZ[X]", RZ, ("X",), (XZ,))

    det_not_one_control = _negative_control(
        "sln-driver-negative-det-not-one",
        qq_sl4_case.id,
        "matrix determinant is not one",
        merge(qq_sl4_case, (;
            matrix = 2 * qq_sl4_case.matrix,
        )),
    )

    unsupported_ring_control = _negative_control(
        "sln-driver-negative-unsupported-coefficient-ring",
        qq_sl4_case.id,
        "coefficient ring is not a field-backed polynomial ring",
        merge(qq_sl4_case, (;
            ring_constructor = zz_ring_constructor,
            ring = zz_ring,
        )),
    )

    corrupt_peel_control = _negative_control(
        "sln-driver-negative-corrupt-peel-expectation",
        gf2_sl4_case.id,
        "expected peel count disagrees with the recorded peel steps",
        merge(gf2_sl4_case, (;
            expected_peel_count = 2,
        )),
    )

    unknown_staged_reason_control = _negative_control(
        "sln-driver-negative-unknown-staged-reason",
        staged_missing_case.id,
        "staged reason code is not part of the supported reason vocabulary",
        merge(staged_missing_case, (;
            staged_reason_codes = (:unknown_reason_code,),
        )),
    )

    false_mainline_control = _negative_control(
        "sln-driver-negative-false-mainline-support",
        staged_missing_case.id,
        "missing final SL3 route evidence is mislabeled as issue 186 mainline support",
        merge(staged_missing_case, (;
            support_role = :issue186_mainline,
            expected_status = :supported,
            staged_reason_codes = (),
        )),
    )

    return (;
        cases = [
            legacy_case,
            gf2_sl4_case,
            gf2_sl5_case,
            qq_sl4_case,
            staged_missing_case,
        ],
        negative_controls = [
            det_not_one_control,
            unsupported_ring_control,
            corrupt_peel_control,
            unknown_staged_reason_control,
            false_mainline_control,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
