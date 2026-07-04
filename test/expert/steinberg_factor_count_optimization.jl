using Test
using Suslin
using Oscar

@testset "Steinberg canonical elementary factor records" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    coefficient = x + y + one(R)

    elementary_factor = elementary_matrix(3, 1, 2, coefficient, R)
    elementary_record = Suslin._canonical_elementary_factor_record(elementary_factor)

    @test elementary_record.kind == :elementary
    @test elementary_record.n == 3
    @test Suslin._same_base_ring(elementary_record.ring, R)
    @test elementary_record.row == 1
    @test elementary_record.col == 2
    @test elementary_record.coefficient == coefficient
    @test Suslin._elementary_factor_record_matrix(elementary_record) == elementary_factor

    zero_elementary_factor = elementary_matrix(3, 1, 2, zero(R), R)
    identity_record = Suslin._canonical_elementary_factor_record(zero_elementary_factor)

    @test identity_record.kind == :identity
    @test identity_record.n == 3
    @test Suslin._same_base_ring(identity_record.ring, R)
    @test Suslin._elementary_factor_record_matrix(identity_record) == zero_elementary_factor

    nonsquare_factor = zero_matrix(R, 2, 3)
    bad_diagonal = identity_matrix(R, 3)
    bad_diagonal[2, 2] = x
    two_offdiagonal = identity_matrix(R, 3)
    two_offdiagonal[1, 2] = x
    two_offdiagonal[2, 3] = y

    @test_throws DimensionMismatch Suslin._canonical_elementary_factor_record(nonsquare_factor)
    @test_throws ArgumentError Suslin._canonical_elementary_factor_record(bad_diagonal)
    @test_throws ArgumentError Suslin._canonical_elementary_factor_record(two_offdiagonal)
end
