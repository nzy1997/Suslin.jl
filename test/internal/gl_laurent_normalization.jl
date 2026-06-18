using Test
using Suslin
using Oscar

@testset "Laurent GL_n determinant classification and normalization" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])

    determinant_one = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])
    one_normalization = normalize_laurent_gl_matrix(determinant_one)
    @test classify_laurent_determinant(determinant_one).classification == :one
    @test one_normalization.determinant_classification == :one
    @test one_normalization.normalized_matrix == determinant_one
    @test one_normalization.correction.kind == :identity
    @test verify_laurent_gl_normalization(determinant_one, one_normalization)

    monomial_unit = matrix(R, [
        x^-1 * y one(R);
        zero(R) one(R)
    ])
    monomial_profile = classify_laurent_determinant(monomial_unit)
    @test monomial_profile.classification == :laurent_monomial_unit
    @test monomial_profile.monomial_exponents == (-1, 1)
    monomial_normalization = normalize_laurent_gl_matrix(monomial_unit)
    @test monomial_normalization.determinant_classification == :laurent_monomial_unit
    @test det(monomial_normalization.normalized_matrix) == one(R)
    @test monomial_normalization.correction.kind == :left_diagonal_determinant_correction
    @test monomial_normalization.correction.side == :left
    @test verify_laurent_gl_normalization(monomial_unit, monomial_normalization)
    @test monomial_normalization.correction.factor * monomial_normalization.normalized_matrix == monomial_unit

    tampered_factor = copy(monomial_normalization.correction.factor)
    tampered_factor[1, 1] = one(R)
    tampered = merge(
        monomial_normalization,
        (;
            correction = merge(
                monomial_normalization.correction,
                (; factor = tampered_factor),
            ),
        ),
    )
    @test !verify_laurent_gl_normalization(monomial_unit, tampered)

    Q, (t,) = suslin_laurent_polynomial_ring(QQ, ["t"])
    transposition = matrix(Q, [
        zero(Q) one(Q);
        one(Q) zero(Q)
    ])
    sign_profile = classify_laurent_determinant(transposition)
    @test sign_profile.classification == :permutation_sign_unit
    sign_normalization = normalize_laurent_gl_matrix(transposition)
    @test sign_normalization.determinant_classification == :permutation_sign_unit
    @test det(sign_normalization.normalized_matrix) == one(Q)
    @test verify_laurent_gl_normalization(transposition, sign_normalization)

    non_unit = matrix(R, [
        x + one(R) zero(R);
        zero(R) one(R)
    ])
    non_unit_profile = classify_laurent_determinant(non_unit)
    @test non_unit_profile.classification == :non_unit
    err = try
        normalize_laurent_gl_matrix(non_unit)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("unsupported Laurent GL_n determinant", sprint(showerror, err))
    @test occursin("outside the staged SL_n factorization path", sprint(showerror, err))

    Z4, _ = residue_ring(ZZ, 4)
    S, (u,) = laurent_polynomial_ring(Z4, ["u"])
    monomial_nonunit = matrix(S, [S(Z4(2)) * u;;])
    monomial_nonunit_profile = classify_laurent_determinant(monomial_nonunit)
    @test monomial_nonunit_profile.classification == :non_unit
    monomial_nonunit_err = try
        normalize_laurent_gl_matrix(monomial_nonunit)
        nothing
    catch caught
        caught
    end
    @test monomial_nonunit_err isa ArgumentError
    @test occursin("unsupported Laurent GL_n determinant", sprint(showerror, monomial_nonunit_err))

    other_unit_value = one(S) + S(Z4(2)) * u
    other_unit = matrix(S, [other_unit_value;;])
    other_unit_profile = classify_laurent_determinant(other_unit)
    @test other_unit_profile.classification == :other_unit
    other_unit_err = try
        normalize_laurent_gl_matrix(other_unit)
        nothing
    catch caught
        caught
    end
    @test other_unit_err isa ArgumentError
    @test occursin("unsupported Laurent GL_n determinant", sprint(showerror, other_unit_err))
    @test occursin("non-monomial units are outside the staged SL_n factorization path", sprint(showerror, other_unit_err))
end
