using Suslin
using Test
using Oscar

const PARK_WOODBURN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

function _captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "public factorization driver shell" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    supported = matrix(R, [
        one(R)      one(R) + X       zero(R);
        X           one(R) + X + X^2 zero(R);
        zero(R)     zero(R)          one(R)
    ])
    factors = elementary_factorization(supported)
    @test verify_factorization(supported, factors)
    supported_cert = Suslin._polynomial_factorization_route_certificate(supported)
    @test supported_cert.route == :fast_local_sl3
    @test Suslin._verify_polynomial_factorization_route_certificate(supported_cert)
    @test factors == supported_cert.factors

    nonfixture_p = X^2 + X + one(R)
    nonfixture_q = X
    nonfixture_local_sl3 = matrix(R, [
        nonfixture_p               nonfixture_q zero(R);
        X + one(R) + nonfixture_p X + one(R)   zero(R);
        zero(R)                    zero(R)      one(R)
    ])
    nonfixture_factors = elementary_factorization(nonfixture_local_sl3)
    @test verify_factorization(nonfixture_local_sl3, nonfixture_factors)
    nonfixture_cert = Suslin._polynomial_factorization_route_certificate(nonfixture_local_sl3)
    @test nonfixture_cert.route == :fast_local_sl3
    @test nonfixture_cert.evidence isa Suslin.SL3LocalRealizationCertificate
    @test nonfixture_cert.evidence.branch == :murthy_q0_nonunit_bezout_resultant
    @test nonfixture_factors == nonfixture_cert.factors

    larger_sl = identity_matrix(R, 4)
    larger_factors = elementary_factorization(larger_sl)
    @test isempty(larger_factors)
    @test verify_factorization(larger_sl, larger_factors)
    larger_cert = Suslin._polynomial_factorization_route_certificate(larger_sl)
    @test larger_cert.route == :disjoint_local_blocks
    @test Suslin._verify_polynomial_factorization_route_certificate(larger_cert)
    @test larger_factors == larger_cert.factors

    unsupported_larger = identity_matrix(R, 4)
    unsupported_larger[1, 4] = X
    larger_err = _captured_error(() -> elementary_factorization(unsupported_larger))
    @test larger_err isa ArgumentError
    @test occursin("staged SL_n to local SL_3 reduction failure", sprint(showerror, larger_err))
    @test !occursin("currently supports only 3x3 matrices", sprint(showerror, larger_err))

    nonsquare = zero_matrix(R, 3, 4)
    nonsquare_err = _captured_error(() -> elementary_factorization(nonsquare))
    @test nonsquare_err isa ArgumentError
    @test occursin("A must be square", sprint(showerror, nonsquare_err))

    undersized = identity_matrix(R, 2)
    undersized_err = _captured_error(() -> elementary_factorization(undersized))
    @test undersized_err isa ArgumentError
    @test occursin("size at least 3", sprint(showerror, undersized_err))

    unsupported_base_ring = identity_matrix(ZZ, 3)
    unsupported_ring_err = _captured_error(() -> elementary_factorization(unsupported_base_ring))
    @test unsupported_ring_err isa ArgumentError
    @test occursin("outside the supported exact polynomial or Laurent polynomial factorization path", sprint(showerror, unsupported_ring_err))

    determinant_not_one = matrix(R, [
        X + one(R) zero(R) zero(R) zero(R);
        zero(R)    one(R) zero(R) zero(R);
        zero(R)    zero(R) one(R) zero(R);
        zero(R)    zero(R) zero(R) one(R)
    ])
    determinant_err = _captured_error(() -> elementary_factorization(determinant_not_one))
    @test determinant_err isa ArgumentError
    @test occursin("determinant/unit precondition", sprint(showerror, determinant_err))
    @test occursin("outside the staged SL_n factorization path", sprint(showerror, determinant_err))
    @test !occursin("SL_n reduction layer", sprint(showerror, determinant_err))
    @test !occursin("currently supports only 3x3 matrices", sprint(showerror, determinant_err))

    S, _ = Oscar.polynomial_ring(QQ, ["X", "Y"])
    multivariate_sl3 = identity_matrix(S, 3)
    multivariate_err = _captured_error(() -> elementary_factorization(multivariate_sl3))
    @test multivariate_err isa ArgumentError
    @test occursin("missing Quillen/local realizability witness", sprint(showerror, multivariate_err))

    nonlocal_sl3 = matrix(R, [
        one(R)  zero(R) X;
        zero(R) one(R)  zero(R);
        zero(R) zero(R) one(R)
    ])
    nonlocal_err = _captured_error(() -> elementary_factorization(nonlocal_sl3))
    @test nonlocal_err isa ArgumentError
    @test occursin("staged reduction to the supported univariate local SL_3 slice", sprint(showerror, nonlocal_err))

    L, (x,) = suslin_laurent_polynomial_ring(GF(2), ["x"])
    normalizable_laurent = matrix(L, [
        x       zero(L) zero(L);
        zero(L) one(L)  zero(L);
        zero(L) zero(L) one(L)
    ])
    laurent_err = _captured_error(() -> elementary_factorization(normalizable_laurent))
    @test laurent_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, laurent_err))
    @test occursin("determinant-correction/driver path cannot yet return elementary factors that reconstruct the original input", sprint(showerror, laurent_err))

    non_normalizable_laurent = matrix(L, [
        x + one(L) zero(L) zero(L);
        zero(L)    one(L)  zero(L);
        zero(L)    zero(L) one(L)
    ])
    non_normalizable_err = _captured_error(() -> elementary_factorization(non_normalizable_laurent))
    @test non_normalizable_err isa ArgumentError
    @test occursin("unsupported Laurent GL_n determinant", sprint(showerror, non_normalizable_err))
    @test occursin("outside the staged SL_n factorization path", sprint(showerror, non_normalizable_err))
    @test !occursin("SL_n reduction layer", sprint(showerror, non_normalizable_err))
    @test !occursin("currently supports only 3x3 matrices", sprint(showerror, non_normalizable_err))

    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_DRIVER_CATALOG_PATH)
    end
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()
    recursive_supported = entries["pw-poly-recursive-column-peel-sl3-qq"].matrix
    recursive_factors = elementary_factorization(recursive_supported)
    @test verify_factorization(recursive_supported, recursive_factors)
    recursive_cert = Suslin._polynomial_factorization_route_certificate(recursive_supported)
    @test recursive_cert.route == :polynomial_column_peel
    @test recursive_factors == recursive_cert.factors
    @test recursive_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test Suslin._verify_polynomial_column_peel_certificate(recursive_cert.evidence)

    quillen_supported = entries["quillen-patched-substitution-witness-qq"].matrix
    quillen_factors = nothing
    quillen_supported_err = _captured_error(() -> begin
        quillen_factors = elementary_factorization(quillen_supported)
        nothing
    end)
    @test quillen_supported_err === nothing
    if quillen_supported_err === nothing
        @test verify_factorization(quillen_supported, quillen_factors)
        quillen_cert = Suslin._polynomial_factorization_route_certificate(quillen_supported)
        @test quillen_cert.route == :quillen_patch
        @test quillen_factors == quillen_cert.factors
        @test quillen_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
        @test Suslin._verify_polynomial_quillen_patch_route_adapter(quillen_cert.evidence)
    end

    S = base_ring(quillen_supported)
    SX, Sr, Sg = collect(gens(S))
    quillen_unsupported = elementary_matrix(
        3,
        1,
        2,
        SX + Sr^2 * Sg + Sg + Sr + one(S),
        S,
    )
    quillen_unsupported_err =
        _captured_error(() -> elementary_factorization(quillen_unsupported))
    @test quillen_unsupported_err isa ArgumentError
    @test occursin("missing Quillen/local realizability witness", sprint(showerror, quillen_unsupported_err))

    recursive_unsupported = entries["pw-poly-recursive-column-peel-gf2"].matrix
    recursive_unsupported_err =
        _captured_error(() -> elementary_factorization(recursive_unsupported))
    @test recursive_unsupported_err isa ArgumentError
    @test occursin(
        "missing Quillen/local realizability witness",
        sprint(showerror, recursive_unsupported_err),
    )
end
