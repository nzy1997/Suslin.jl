module ParkWoodburnMainlineAcceptanceFixtureCatalog

using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_sl3_driver_cases.jl")
const ECP_MAINLINE_CATALOG_PATH =
    joinpath(@__DIR__, "ecp_mainline_cases.jl")
const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_sln_driver_cases.jl")

const PARK_WOODBURN_SECTION_2_REF = "refs/arXiv-alg-geom9405003v1 Section 2"
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"

if !isdefined(@__MODULE__, :ParkWoodburnSL3DriverFixtureCatalog)
    include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
end
if !isdefined(@__MODULE__, :ECPMainlineFixtureCatalog)
    include(ECP_MAINLINE_CATALOG_PATH)
end
if !isdefined(@__MODULE__, :ParkWoodburnSLnDriverFixtureCatalog)
    include(PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH)
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

function _public_route(route_marker::Symbol, status::Symbol)
    return (;
        entrypoint = :elementary_factorization,
        issue_id = "#187",
        route_marker = route_marker,
        status = status,
    )
end

function _parent_issue_coverage(;
    issue181,
    issue182,
    issue183,
    issue184,
    issue185,
    issue186,
)
    return (; issue181, issue182, issue183, issue184, issue185, issue186)
end

function _acceptance_metadata(boundary::Symbol, parent_issue_coverage, note::AbstractString)
    return (;
        boundary = boundary,
        parent_issue_coverage = parent_issue_coverage,
        note = note,
    )
end

function _case(;
    id,
    entry_class,
    expected_status,
    public_route,
    ring_constructor,
    ring,
    matrix,
    determinant_metadata,
    source_refs,
    upstream_issue_ids,
    upstream_evidence,
    acceptance_metadata,
    missing_evidence = nothing,
    staged_reason = nothing,
)
    entry = (;
        id = id,
        entry_class = entry_class,
        expected_status = expected_status,
        public_route = public_route,
        ring_constructor = ring_constructor,
        ring = ring,
        matrix = matrix,
        determinant_metadata = determinant_metadata,
        source_refs = source_refs,
        upstream_issue_ids = upstream_issue_ids,
        upstream_evidence = upstream_evidence,
        acceptance_metadata = acceptance_metadata,
    )
    missing_evidence === nothing || (entry = merge(entry, (; missing_evidence = missing_evidence)))
    staged_reason === nothing || (entry = merge(entry, (; staged_reason = staged_reason)))
    return entry
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (; id = id, base_case_id = base_case_id, reason = reason))
end

function _determinant_metadata(matrix_value; certified = true)
    R = base_ring(matrix_value)
    return (;
        expected = :one,
        certified = certified,
        value = det(matrix_value),
        unit = one(R),
    )
end

function _sl3_cases_by_id()
    return Base.invokelatest(getfield(ParkWoodburnSL3DriverFixtureCatalog, :cases_by_id))
end

function _ecp_cases_by_id()
    return Base.invokelatest(getfield(ECPMainlineFixtureCatalog, :cases_by_id))
end

function _sln_cases_by_id()
    return Base.invokelatest(getfield(ParkWoodburnSLnDriverFixtureCatalog, :cases_by_id))
end

function catalog()
    sl3_cases = _sl3_cases_by_id()
    ecp_cases = _ecp_cases_by_id()
    sln_cases = _sln_cases_by_id()

    sl3_multivariate_base = sl3_cases["sl3-driver-multivariate-monic-special-form-qq"]
    readme_base = sl3_cases["sl3-driver-univariate-fast-local-qq"]
    sln_recursive_base = sln_cases["sln-driver-sl4-gf2-ecp-mainline"]
    staged_base = sln_cases["sln-driver-staged-missing-final-sl3-qq"]
    ecp_mainline_base = ecp_cases["ecp-mainline-gf2-hard-slice"]
    final_sl3_base = sl3_cases["sl3-driver-quillen-mainline-evidence-gf2"]

    sl3_multivariate_case = _case(
        id = "pw-mainline-sl3-multivariate-issue184-qq",
        entry_class = :issue184_sl3_multivariate,
        expected_status = :mainline_accepted,
        public_route = _public_route(:issue184_evidence_backed_sl3, :mainline_accepted),
        ring_constructor = sl3_multivariate_base.ring_constructor,
        ring = sl3_multivariate_base.ring,
        matrix = sl3_multivariate_base.matrix,
        determinant_metadata = _determinant_metadata(sl3_multivariate_base.matrix),
        source_refs = (
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 270 final acceptance reuses the issue 184 multivariate SL_3 driver witness",
        ),
        upstream_issue_ids = ("#184",),
        upstream_evidence = (; driver_case_id = sl3_multivariate_base.id),
        acceptance_metadata = _acceptance_metadata(
            :supported_public_example,
            _parent_issue_coverage(
                issue181 = :covered,
                issue182 = :covered,
                issue183 = :covered,
                issue184 = :covered,
                issue185 = :not_applicable,
                issue186 = :not_applicable,
            ),
            "Evidence-backed multivariate SL_3 acceptance shell for the final ordinary-polynomial mainline claim",
        ),
    )

    readme_style_case = _case(
        id = "pw-mainline-readme-ordinary-polynomial-qq",
        entry_class = :readme_public_example,
        expected_status = :mainline_accepted,
        public_route = _public_route(:issue184_readme_public_example, :mainline_accepted),
        ring_constructor = readme_base.ring_constructor,
        ring = readme_base.ring,
        matrix = readme_base.matrix,
        determinant_metadata = _determinant_metadata(readme_base.matrix),
        source_refs = (
            PARK_WOODBURN_SECTION_2_REF,
            PARK_WOODBURN_SECTION_5_REF,
            "README Example ordinary-polynomial public shell",
        ),
        upstream_issue_ids = ("#184",),
        upstream_evidence = (; driver_case_id = readme_base.id),
        acceptance_metadata = _acceptance_metadata(
            :supported_readme_example,
            _parent_issue_coverage(
                issue181 = :covered,
                issue182 = :covered,
                issue183 = :covered,
                issue184 = :covered,
                issue185 = :not_applicable,
                issue186 = :not_applicable,
            ),
            "README-style determinant-one example counted as a public ordinary-polynomial acceptance shell",
        ),
    )

    sln_recursive_case = _case(
        id = "pw-mainline-sln-recursive-issue185-186-gf2",
        entry_class = :issue185_186_sln_recursive,
        expected_status = :mainline_accepted,
        public_route = _public_route(:issue186_recursive_mainline, :mainline_accepted),
        ring_constructor = sln_recursive_base.ring_constructor,
        ring = sln_recursive_base.ring,
        matrix = sln_recursive_base.matrix,
        determinant_metadata = _determinant_metadata(sln_recursive_base.matrix),
        source_refs = (
            PARK_WOODBURN_SECTION_3_REF,
            PARK_WOODBURN_SECTION_4_REF,
            PARK_WOODBURN_SECTION_5_REF,
            "Issue 270 final acceptance reuses issue 185 ECP and issue 186 recursive SL_n evidence",
        ),
        upstream_issue_ids = ("#184", "#185", "#186"),
        upstream_evidence = (;
            sln_case_id = sln_recursive_base.id,
            ecp_case_id = ecp_mainline_base.id,
            final_sl3_case_id = final_sl3_base.id,
        ),
        acceptance_metadata = _acceptance_metadata(
            :supported_recursive_mainline,
            _parent_issue_coverage(
                issue181 = :covered,
                issue182 = :covered,
                issue183 = :covered,
                issue184 = :covered,
                issue185 = :covered,
                issue186 = :covered,
            ),
            "Recursive ordinary-polynomial SL_n acceptance shell with replayed ECP peel steps and final issue 184 SL_3 evidence",
        ),
    )

    staged_missing_evidence_case = _case(
        id = "pw-mainline-staged-missing-evidence-qq",
        entry_class = :staged_missing_evidence_boundary,
        expected_status = :staged,
        public_route = _public_route(:staged_missing_issue184_final_route, :staged),
        ring_constructor = staged_base.ring_constructor,
        ring = staged_base.ring,
        matrix = staged_base.matrix,
        determinant_metadata = _determinant_metadata(staged_base.matrix),
        source_refs = (
            PARK_WOODBURN_SECTION_2_REF,
            PARK_WOODBURN_SECTION_4_REF,
            "Issue 270 staged determinant-one recursive boundary without final issue 184 SL_3 evidence",
        ),
        upstream_issue_ids = ("#186",),
        upstream_evidence = (; sln_case_id = staged_base.id),
        acceptance_metadata = _acceptance_metadata(
            :staged_missing_evidence,
            _parent_issue_coverage(
                issue181 = :not_applicable,
                issue182 = :not_applicable,
                issue183 = :not_applicable,
                issue184 = :staged,
                issue185 = :covered,
                issue186 = :covered,
            ),
            "Determinant-one recursive SL_n acceptance boundary that remains staged until final issue 184 evidence is present",
        ),
        missing_evidence = (:final_sl3_case_id,),
        staged_reason = "missing verified final issue 184 SL_3 evidence for the #187 ordinary-polynomial acceptance claim",
    )

    R_readme = readme_style_case.ring.object
    X_readme = readme_style_case.ring.generators[1]
    det_not_one_matrix = matrix(R_readme, [
        one(R_readme) + X_readme  one(R_readme) + X_readme zero(R_readme);
        X_readme                  one(R_readme) + X_readme^2 zero(R_readme);
        zero(R_readme)            zero(R_readme)            one(R_readme)
    ])
    det_not_one_control = _negative_control(
        "pw-mainline-negative-det-not-one",
        readme_style_case.id,
        "supported mainline acceptance must not record a determinant-not-one matrix",
        merge(readme_style_case, (;
            matrix = det_not_one_matrix,
            determinant_metadata = _determinant_metadata(det_not_one_matrix; certified = false),
        )),
    )

    RZ, (XZ,) = Oscar.polynomial_ring(ZZ, ["X"])
    zz_ring = _ring_metadata("ZZ[X]", RZ, ("X",), (XZ,))
    zz_constructor = _ordinary_ring_constructor("ZZ", ("X",))
    unsupported_ring_matrix = matrix(RZ, [
        one(RZ)      one(RZ) + XZ      zero(RZ);
        XZ           one(RZ) + XZ + XZ^2 zero(RZ);
        zero(RZ)     zero(RZ)         one(RZ)
    ])
    unsupported_ring_control = _negative_control(
        "pw-mainline-negative-unsupported-coefficient-ring",
        readme_style_case.id,
        "supported mainline acceptance must stay on exact field-backed ordinary polynomial rings",
        merge(readme_style_case, (;
            ring_constructor = zz_constructor,
            ring = zz_ring,
            matrix = unsupported_ring_matrix,
            determinant_metadata = _determinant_metadata(unsupported_ring_matrix),
        )),
    )

    missing_evidence_control = _negative_control(
        "pw-mainline-negative-missing-evidence",
        sln_recursive_case.id,
        "mainline acceptance cannot claim recursive support while omitting the final issue 184 evidence id",
        merge(sln_recursive_case, (;
            upstream_evidence = (;
                sln_case_id = sln_recursive_base.id,
                ecp_case_id = ecp_mainline_base.id,
            ),
        )),
    )

    return (;
        cases = [
            sl3_multivariate_case,
            sln_recursive_case,
            readme_style_case,
            staged_missing_evidence_case,
        ],
        negative_controls = [
            det_not_one_control,
            unsupported_ring_control,
            missing_evidence_control,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
