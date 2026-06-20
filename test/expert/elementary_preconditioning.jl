using Test
using Suslin
using Oscar

@testset "side-aware elementary preconditioning" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    A = matrix(R, [
        one(R)      x                 y^-1;
        y           one(R) + x * y    x^-1;
        zero(R)     y^-1              one(R) + x
    ])

    left_step = elementary_preconditioning_step(A, :left, 2, 1, x^-1)
    expected_left_factor = elementary_matrix(3, 2, 1, x^-1, R)
    expected_left_matrix = expected_left_factor * A

    @test left_step.side == :left
    @test left_step.target == 2
    @test left_step.source == 1
    @test left_step.coefficient == x^-1
    @test left_step.factor == expected_left_factor
    @test left_step.transformed_matrix == expected_left_matrix

    right_step = elementary_preconditioning_step(left_step.transformed_matrix, :right, 3, 1, y)
    expected_right_factor = elementary_matrix(3, 1, 3, y, R)
    expected_right_matrix = left_step.transformed_matrix * expected_right_factor

    @test right_step.side == :right
    @test right_step.target == 3
    @test right_step.source == 1
    @test right_step.coefficient == y
    @test right_step.factor == expected_right_factor
    @test right_step.transformed_matrix == expected_right_matrix

    final_step = elementary_preconditioning_step(right_step.transformed_matrix, :left, 1, 3, x * y^-1)
    expected_final_factor = elementary_matrix(3, 1, 3, x * y^-1, R)
    expected_final_matrix = expected_final_factor * right_step.transformed_matrix
    steps = [left_step, right_step, final_step]

    @test final_step.factor == expected_final_factor
    @test final_step.transformed_matrix == expected_final_matrix
    @test replay_elementary_preconditioning(A, steps) == expected_final_matrix
    @test verify_elementary_preconditioning(A, steps, expected_final_matrix)

    swapped_steps = [
        (; step..., side = step.side == :left ? :right : :left)
        for step in steps
    ]
    @test !verify_elementary_preconditioning(A, swapped_steps, expected_final_matrix)
end

@testset "elementary preconditioning validation" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    A = matrix(R, [
        one(R)  x        zero(R);
        y       one(R)   x^-1;
        zero(R) y^-1     one(R)
    ])

    @test_throws ArgumentError elementary_preconditioning_step(A, :middle, 1, 2, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :left, 1, 1, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :left, 0, 1, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :right, 4, 1, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :right, 1, 2, u)

    wrong_ring_factor = elementary_matrix(3, 1, 2, u, S)
    wrong_size_factor = elementary_matrix(2, 1, 2, one(R), R)

    @test_throws ArgumentError replay_elementary_preconditioning(A, [(; side = :left, factor = wrong_ring_factor)])
    @test_throws ArgumentError replay_elementary_preconditioning(A, [(; side = :right, factor = wrong_size_factor)])
    @test_throws ArgumentError replay_elementary_preconditioning(A, [(; factor = identity_matrix(R, 3))])

    @test !verify_elementary_preconditioning(A, [(; side = :left, factor = wrong_ring_factor)], A)
    @test !verify_elementary_preconditioning(A, [(; side = :right, factor = wrong_size_factor)], A)
end
