using Suslin
using Test
using Oscar

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

    larger_sl = identity_matrix(R, 4)
    larger_err = _captured_error(() -> elementary_factorization(larger_sl))
    @test larger_err isa ArgumentError
    @test occursin("SL_n reduction layer", sprint(showerror, larger_err))
    @test occursin("not yet implemented", sprint(showerror, larger_err))
    @test !occursin("currently supports only 3x3 matrices", sprint(showerror, larger_err))

    nonsquare = zero_matrix(R, 3, 4)
    nonsquare_err = _captured_error(() -> elementary_factorization(nonsquare))
    @test nonsquare_err isa ArgumentError
    @test occursin("A must be square", sprint(showerror, nonsquare_err))

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

    L, (x,) = suslin_laurent_polynomial_ring(GF(2), ["x"])
    normalizable_laurent = matrix(L, [
        x       zero(L) zero(L);
        zero(L) one(L)  zero(L);
        zero(L) zero(L) one(L)
    ])
    laurent_err = _captured_error(() -> elementary_factorization(normalizable_laurent))
    @test laurent_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, laurent_err))
    @test occursin("Laurent SL_n reduction layer", sprint(showerror, laurent_err))
    @test occursin("not yet implemented", sprint(showerror, laurent_err))

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
end
