using Suslin
using Test
using Oscar

const PARK_WOODBURN_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")
const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")
const PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_mainline_acceptance_cases.jl")

function _pw_acceptance_result_or_error(A)
    factors = nothing
    try
        factors = elementary_factorization(A)
        return factors, nothing
    catch err
        return factors, err
    end
end

function _pw_assert_public_issue186_recursive_acceptance(A, expected_step_count::Int)
    factors, err = _pw_acceptance_result_or_error(A)
    @test err === nothing
    @test factors !== nothing
    @test verify_factorization(A, factors)
    cert = Suslin._polynomial_factorization_route_certificate(A)
    @test cert.route == :polynomial_column_peel
    @test cert.status == :supported
    @test factors == cert.factors
    @test cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test length(cert.evidence.peel_steps) == expected_step_count
    @test cert.evidence.descent_metadata.descent_dimensions ==
          tuple((nrows(A):-1:3)...)
    @test cert.evidence.mainline_support_metadata.issue_id == "#186"
    @test cert.evidence.mainline_support_metadata.marker == :issue186_mainline
    @test cert.evidence.mainline_support_metadata.supported
    @test cert.evidence.mainline_support_metadata.peel_steps_ecp_verified
    @test cert.evidence.mainline_support_metadata.final_route_issue184_ok
    @test cert.evidence.final_route_provenance ==
          :issue184_evidence_backed_sl3
    @test all(step -> Suslin.verify_ecp_column_reduction(step.ecp_evidence),
        cert.evidence.peel_steps)
    @test cert.evidence.final_certificate.route == :quillen_patch
    @test Suslin._verify_polynomial_column_peel_certificate(cert.evidence)
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)

    corrupted = copy(factors)
    corrupted[1] =
        corrupted[1] *
        elementary_matrix(nrows(A), 1, 2, one(base_ring(A)), base_ring(A))
    @test !verify_factorization(A, corrupted)

    return cert
end

function _pw_assert_issue184_sl3_public_acceptance(entry)
    A = entry.matrix
    @test nrows(A) == 3
    @test length(collect(gens(base_ring(A)))) > 1
    factors, err = _pw_acceptance_result_or_error(A)
    @test err === nothing
    @test factors !== nothing
    @test verify_factorization(A, factors)
    cert = Suslin._polynomial_factorization_route_certificate(A)
    @test cert.route == :quillen_patch
    @test cert.status == :supported
    @test factors == cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)
    @test cert.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence ||
          cert.evidence isa Suslin.PolynomialSL3SuppliedQuillenRouteEvidence ||
          cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
    if cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
        issue184_cert = Suslin._polynomial_sl3_supplied_quillen_route_certificate(A)
        @test issue184_cert.route == :quillen_patch
        @test issue184_cert.status == :supported
        @test Suslin._verify_polynomial_factorization_route_certificate(issue184_cert)
        @test issue184_cert.evidence isa Suslin.PolynomialSL3SuppliedQuillenRouteEvidence
        @test issue184_cert.evidence.replay_metadata.driver_issue_id == "#184"
        @test issue184_cert.factors == factors
        @test issue184_cert.evidence.quillen_route_adapter.target == cert.evidence.target
        @test issue184_cert.evidence.quillen_route_adapter.target_matrix ==
              cert.evidence.target_matrix
        @test issue184_cert.evidence.quillen_route_adapter.product == cert.evidence.product
        @test issue184_cert.evidence.quillen_route_adapter.global_elementary_factors ==
              cert.evidence.global_elementary_factors
    end
    @test entry.public_route.issue_id == "#187"
    @test "#184" in entry.upstream_issue_ids
    return cert
end

function _pw_assert_readme_public_acceptance(entry)
    A = entry.matrix
    factors, err = _pw_acceptance_result_or_error(A)
    @test err === nothing
    @test factors !== nothing
    @test verify_factorization(A, factors)
    cert = Suslin._polynomial_factorization_route_certificate(A)
    @test cert.status == :supported
    @test factors == cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)
    @test entry.entry_class == :readme_public_example
    @test entry.public_route.issue_id == "#187"
    return cert
end

function _pw_assert_issue187_recursive_catalog_acceptance(entry, expected_step_count::Int)
    cert = _pw_assert_public_issue186_recursive_acceptance(entry.matrix, expected_step_count)
    @test entry.entry_class == :issue185_186_sln_recursive
    @test entry.public_route.route_marker == :issue186_recursive_mainline
    @test entry.public_route.issue_id == "#187"
    @test "#184" in entry.upstream_issue_ids
    @test "#185" in entry.upstream_issue_ids
    @test "#186" in entry.upstream_issue_ids
    @test hasproperty(entry.upstream_evidence, :ecp_case_id)
    @test hasproperty(entry.upstream_evidence, :final_sl3_case_id)
    return cert
end

@testset "public Park-Woodburn polynomial factorization acceptance" begin
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_ACCEPTANCE_CATALOG_PATH)
    end
    if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        include(PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH)
    end
    if !isdefined(Main, :ParkWoodburnMainlineAcceptanceFixtureCatalog)
        include(PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH)
    end

    catalog = ParkWoodburnPolynomialFixtureCatalog.catalog()
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()
    sln_entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
    mainline_entries = ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()

    sl3_mainline =
        mainline_entries["pw-mainline-sl3-multivariate-issue184-qq"]
    _pw_assert_issue184_sl3_public_acceptance(sl3_mainline)

    recursive_mainline =
        mainline_entries["pw-mainline-sln-recursive-issue185-186-gf2"]
    _pw_assert_issue187_recursive_catalog_acceptance(recursive_mainline, 1)

    readme_mainline =
        mainline_entries["pw-mainline-readme-ordinary-polynomial-qq"]
    _pw_assert_readme_public_acceptance(readme_mainline)

    fast_local = entries["pw-poly-univariate-sl3-fast-local-qq"].matrix
    fast_factors, fast_err = _pw_acceptance_result_or_error(fast_local)
    @test fast_err === nothing
    @test fast_factors !== nothing
    @test verify_factorization(fast_local, fast_factors)
    fast_cert = Suslin._polynomial_factorization_route_certificate(fast_local)
    @test fast_cert.route == :fast_local_sl3
    @test fast_factors == fast_cert.factors

    recursive = sln_entries["sln-driver-sl4-gf2-ecp-mainline"].matrix
    _pw_assert_public_issue186_recursive_acceptance(recursive, 1)

    recursive_two_step = sln_entries["sln-driver-sl5-gf2-two-step"].matrix
    _pw_assert_public_issue186_recursive_acceptance(recursive_two_step, 2)

    quillen = entries["quillen-patched-substitution-witness-qq"].matrix
    quillen_factors, quillen_err = _pw_acceptance_result_or_error(quillen)
    @test quillen_err === nothing
    @test quillen_factors !== nothing
    @test verify_factorization(quillen, quillen_factors)
    quillen_cert = Suslin._polynomial_factorization_route_certificate(quillen)
    @test quillen_cert.route == :quillen_patch
    @test quillen_factors == quillen_cert.factors
    @test quillen_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test quillen_cert.evidence.quillen_patch isa Suslin.QuillenSuppliedEvidencePatchAssembly
    @test Suslin._verify_polynomial_quillen_patch_route_adapter(quillen_cert.evidence)
    @test Suslin._verify_polynomial_factorization_route_certificate(quillen_cert)

    issue238_R, (issue238_X, issue238_r, issue238_g) =
        Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    issue238_p = issue238_X + issue238_r * issue238_g + one(issue238_R)
    issue238_A = matrix(issue238_R, [
        issue238_p one(issue238_R) zero(issue238_R);
        issue238_X + issue238_r * issue238_g one(issue238_R) zero(issue238_R);
        zero(issue238_R) zero(issue238_R) one(issue238_R)
    ])
    issue238_factors, issue238_err = _pw_acceptance_result_or_error(issue238_A)
    @test issue238_err === nothing
    @test issue238_factors !== nothing
    @test verify_factorization(issue238_A, issue238_factors)
    issue238_cert = Suslin._polynomial_factorization_route_certificate(issue238_A)
    @test issue238_cert.route == :quillen_patch
    @test issue238_cert.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence
    @test issue238_factors == issue238_cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(issue238_cert)
    @test issue238_cert.status == :supported
    issue238_evidence = issue238_cert.evidence
    @test issue238_evidence.context.catalog_metadata.context_issue_id == "#235"
    @test issue238_evidence.context.catalog_metadata.driver_issue_id == "#184"
    @test Suslin._verify_sl3_realization_input_context(issue238_evidence.context)
    @test issue238_evidence.witness_selection.local_form_witness.witness_issue_id == "#236"
    @test Suslin._verify_sl3_local_form_witness_selection(
        issue238_evidence.witness_selection,
    )
    @test issue238_evidence.local_evidence_provider.metadata.provider_issue_id == "#237"
    @test Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        issue238_evidence.local_evidence_provider,
    )
    @test Suslin.verify_quillen_murthy_adapter_consumption(
        issue238_evidence.quillen_consumption,
    )
    @test Suslin._verify_polynomial_quillen_patch_route_adapter(
        issue238_evidence.quillen_route_adapter,
    )
    issue238_route_metadata =
        issue238_evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata
    @test issue238_route_metadata.context_issue_id == "#235"
    @test issue238_route_metadata.witness_issue_id == "#236"
    @test issue238_route_metadata.provider_issue_id == "#237"
    @test issue238_route_metadata.patch_issue_id == "#220"

    corrupted_issue238_factors = copy(issue238_factors)
    corrupted_issue238_factors[1] =
        corrupted_issue238_factors[1] *
        elementary_matrix(3, 1, 2, one(issue238_R), issue238_R)
    @test !verify_factorization(issue238_A, corrupted_issue238_factors)

    missing_issue236_witness =
        elementary_matrix(3, 1, 3, issue238_X, issue238_R) *
        elementary_matrix(3, 2, 1, issue238_r, issue238_R)
    missing_issue236_factors, missing_issue236_err =
        _pw_acceptance_result_or_error(missing_issue236_witness)
    @test missing_issue236_factors === nothing
    @test missing_issue236_err isa ArgumentError
    @test occursin(
        "evidence-backed SL_3 polynomial route",
        sprint(showerror, missing_issue236_err),
    )
    @test occursin("#236 local-form witness", sprint(showerror, missing_issue236_err))

    S = base_ring(quillen)
    X, r, g = collect(gens(S))
    nonfixture_quillen = elementary_matrix(
        3,
        1,
        3,
        X * r + g + one(S),
        S,
    )
    nonfixture_factors, nonfixture_err = _pw_acceptance_result_or_error(nonfixture_quillen)
    @test nonfixture_err === nothing
    @test nonfixture_factors !== nothing
    @test verify_factorization(nonfixture_quillen, nonfixture_factors)
    nonfixture_cert = Suslin._polynomial_factorization_route_certificate(nonfixture_quillen)
    @test nonfixture_cert.route == :quillen_patch
    @test nonfixture_factors == nonfixture_cert.factors
    @test nonfixture_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test nonfixture_cert.evidence.quillen_patch isa Suslin.QuillenSuppliedEvidencePatchAssembly
    @test Suslin._verify_polynomial_factorization_route_certificate(nonfixture_cert)

    negative_entries = Dict(entry.id => entry for entry in catalog.negative_controls)
    det_factors, det_err =
        _pw_acceptance_result_or_error(negative_entries["pw-poly-det-not-one-control"].matrix)
    @test det_factors === nothing
    @test det_err isa ArgumentError
    @test occursin("determinant/unit precondition", sprint(showerror, det_err))
    @test occursin("outside the staged SL_n factorization path", sprint(showerror, det_err))

    outside_factors, outside_err = _pw_acceptance_result_or_error(
        negative_entries["pw-poly-det-one-outside-witness-control"].matrix,
    )
    @test outside_factors === nothing
    @test outside_err isa ArgumentError
    @test occursin(
        "staged reduction to the supported univariate local SL_3 slice",
        sprint(showerror, outside_err),
    )
end
