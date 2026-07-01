using Suslin
using Test
using Oscar

const PARK_WOODBURN_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

function _pw_acceptance_result_or_error(A)
    factors = nothing
    try
        factors = elementary_factorization(A)
        return factors, nothing
    catch err
        return factors, err
    end
end

@testset "public Park-Woodburn polynomial factorization acceptance" begin
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_ACCEPTANCE_CATALOG_PATH)
    end

    catalog = ParkWoodburnPolynomialFixtureCatalog.catalog()
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()

    fast_local = entries["pw-poly-univariate-sl3-fast-local-qq"].matrix
    fast_factors, fast_err = _pw_acceptance_result_or_error(fast_local)
    @test fast_err === nothing
    @test fast_factors !== nothing
    @test verify_factorization(fast_local, fast_factors)
    fast_cert = Suslin._polynomial_factorization_route_certificate(fast_local)
    @test fast_cert.route == :fast_local_sl3
    @test fast_factors == fast_cert.factors

    recursive = entries["pw-poly-recursive-column-peel-sln-block-qq"].matrix
    recursive_factors, recursive_err = _pw_acceptance_result_or_error(recursive)
    @test recursive_err === nothing
    @test recursive_factors !== nothing
    @test verify_factorization(recursive, recursive_factors)
    recursive_cert = Suslin._polynomial_factorization_route_certificate(recursive)
    @test recursive_cert.route == :polynomial_column_peel
    @test recursive_factors == recursive_cert.factors
    @test recursive_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test recursive_cert.evidence.final_certificate.route == :disjoint_local_blocks
    @test Suslin._verify_polynomial_column_peel_certificate(recursive_cert.evidence)
    @test Suslin._verify_polynomial_factorization_route_certificate(recursive_cert)

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
