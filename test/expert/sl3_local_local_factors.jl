using Test
using Suslin
using Oscar

function _local_factor_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _local_factor_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _local_factor_witness(unit, selected_variable, generator)
    R = parent(unit)
    return (;
        context = (;
            kind = :localization_at_maximal_ideal,
            selected_variable,
            maximal_ideal_generators = (generator,),
        ),
        unit,
        residue_unit = one(R),
        residue_inverse = one(R),
        maximal_ideal_generators = (generator,),
        residue_difference_coefficients = (one(R),),
        global_unit = is_unit(unit),
    )
end

@testset "SL3 local elementary factor replay" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    coefficient = X + one(R)
    ordinary_factor = Suslin.elementary_matrix(3, 1, 2, coefficient, R)
    local_factor = Suslin.sl3_local_elementary_factor(1, 2, coefficient, one(R), X)

    @test Suslin.verify_sl3_local_elementary_factor(local_factor)
    @test Suslin.sl3_local_materialize_elementary_factor(local_factor) == ordinary_factor

    ordinary_replay = Suslin.sl3_local_elementary_factor_replay(
        ordinary_factor,
        [local_factor],
        X,
    )
    @test ordinary_replay.mode == :ordinary
    @test ordinary_replay.denominator_product == one(R)
    @test ordinary_replay.materialized_factors == [ordinary_factor]
    @test ordinary_replay.cleared_product == ordinary_factor
    @test Suslin.verify_sl3_local_elementary_factor_replay(ordinary_replay)

    adapted = Suslin.sl3_local_denominator_one_records_from_matrices([ordinary_factor], X)
    @test length(adapted) == 1
    @test Suslin.sl3_local_materialize_elementary_factor(first(adapted)) == ordinary_factor
    adapted_replay = Suslin.sl3_local_elementary_factor_replay(ordinary_factor, adapted, X)
    @test adapted_replay.mode == :ordinary
    @test Suslin.verify_sl3_local_elementary_factor_replay(adapted_replay)

    RU, (u, UX) = Oscar.polynomial_ring(QQ, ["u", "X"])
    denominator = one(RU) + u
    local_unit_witness = _local_factor_witness(denominator, UX, u)
    first_local = Suslin.sl3_local_elementary_factor(
        1,
        2,
        one(RU),
        denominator,
        UX;
        local_unit_witness,
    )
    second_local = Suslin.sl3_local_elementary_factor(
        1,
        2,
        -one(RU),
        denominator,
        UX;
        local_unit_witness,
    )
    identity_target = identity_matrix(RU, 3)
    cleared_replay = Suslin.sl3_local_elementary_factor_replay(
        identity_target,
        [first_local, second_local],
        UX,
    )
    @test cleared_replay.mode == :denominator_cleared
    @test cleared_replay.materialized_factors === nothing
    @test cleared_replay.denominator_product == denominator^2
    @test cleared_replay.cleared_product == denominator^2 * identity_target
    @test Suslin.verify_sl3_local_elementary_factor_replay(cleared_replay)

    nonmaterialized_err =
        _local_factor_error(() -> Suslin.sl3_local_materialize_elementary_factor(first_local))
    @test nonmaterialized_err isa ArgumentError
    @test occursin("cannot materialize", sprint(showerror, nonmaterialized_err))

    corrupted_numerator = Suslin.SL3LocalElementaryFactor(
        first_local.R,
        first_local.n,
        first_local.row,
        first_local.col,
        first_local.numerator + one(RU),
        first_local.denominator,
        first_local.selected_variable,
        first_local.local_unit_witness,
    )
    @test Suslin.verify_sl3_local_elementary_factor(corrupted_numerator)
    bad_numerator_replay = Suslin.SL3LocalElementaryFactorReplay(
        cleared_replay.target,
        [corrupted_numerator, second_local],
        cleared_replay.selected_variable,
        cleared_replay.mode,
        cleared_replay.denominator_product,
        cleared_replay.cleared_product,
        cleared_replay.materialized_factors,
    )
    @test !Suslin.verify_sl3_local_elementary_factor_replay(bad_numerator_replay)

    bad_row = Suslin.SL3LocalElementaryFactor(
        first_local.R,
        first_local.n,
        4,
        first_local.col,
        first_local.numerator,
        first_local.denominator,
        first_local.selected_variable,
        first_local.local_unit_witness,
    )
    @test !Suslin.verify_sl3_local_elementary_factor(bad_row)
    @test_throws ArgumentError Suslin.sl3_local_elementary_factor(
        1,
        1,
        one(RU),
        denominator,
        UX;
        local_unit_witness,
    )
    @test_throws ArgumentError Suslin.sl3_local_elementary_factor(
        1,
        2,
        one(R),
        one(R),
        X + one(R),
    )

    corrupted_witness = merge(
        local_unit_witness,
        (; residue_difference_coefficients = (zero(RU),)),
    )
    @test_throws ArgumentError Suslin.sl3_local_elementary_factor(
        1,
        2,
        one(RU),
        denominator,
        UX;
        local_unit_witness = corrupted_witness,
    )
    bad_witness_record = Suslin.SL3LocalElementaryFactor(
        first_local.R,
        first_local.n,
        first_local.row,
        first_local.col,
        first_local.numerator,
        first_local.denominator,
        first_local.selected_variable,
        corrupted_witness,
    )
    @test !Suslin.verify_sl3_local_elementary_factor(bad_witness_record)

    wrong_denominator_err = _local_factor_error(() -> Suslin.sl3_local_elementary_factor(
        1,
        2,
        one(RU),
        denominator + u,
        UX;
        local_unit_witness,
    ))
    @test wrong_denominator_err isa ArgumentError
    @test occursin("unit does not match expected value", sprint(showerror, wrong_denominator_err))
end
