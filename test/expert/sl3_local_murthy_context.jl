using Test
using Suslin
using Oscar

const SL3_CONTEXT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _context_as_namedtuple(context)
    names = propertynames(context)
    return NamedTuple{names}(Tuple(getproperty(context, name) for name in names))
end

function _captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _local_unit_witness_like(template, unit, coefficients)
    R = parent(unit)
    return merge(
        template,
        (;
            unit,
            residue_unit = one(R),
            residue_inverse = one(R),
            residue_difference_coefficients = coefficients,
            global_unit = is_unit(unit),
        ),
    )
end

@testset "Murthy local input context" begin
    include(SL3_CONTEXT_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    ordinary = by_id["mg-q0-unit-recursion"]
    ordinary_context = Suslin.sl3_local_murthy_input_context(
        ordinary.target,
        ordinary.variable,
    )
    @test ordinary_context.R === ordinary.ring.object
    @test ordinary_context.X == ordinary.variable
    @test ordinary_context.entries == ordinary.entries
    @test ordinary_context.target == ordinary.target
    @test ordinary_context.determinant == one(ordinary.ring.object)
    @test ordinary_context.degree_p == 1
    @test ordinary_context.degree_q == 0
    @test ordinary_context.p0 == one(ordinary.ring.object)
    @test ordinary_context.q0 == one(ordinary.ring.object)
    @test ordinary_context.p_monic == true
    @test ordinary_context.global_units.q0 == true
    @test ordinary_context.local_units.q0 == true
    @test Suslin.verify_sl3_local_murthy_input_context(ordinary_context)
    @test !Suslin.verify_sl3_local_murthy_input_context((; target = nothing))

    unrelated_split_source = by_id["mg-split-lemma-x-square"]
    unrelated_split_witness = first(unrelated_split_source.witnesses)
    @test Suslin._sl3_local_murthy_verify_split_witness(
        unrelated_split_source.ring.object,
        unrelated_split_witness,
    )
    split_nested_context = Suslin.sl3_local_murthy_input_context(
        unrelated_split_source.target,
        unrelated_split_source.variable;
        witness = (; split = unrelated_split_witness),
    )
    @test split_nested_context.split_witness == unrelated_split_witness
    @test Suslin.verify_sl3_local_murthy_input_context(split_nested_context)
    split_direct_context = Suslin.sl3_local_murthy_input_context(
        unrelated_split_source.target,
        unrelated_split_source.variable;
        witness = unrelated_split_witness,
    )
    @test split_direct_context.split_witness == unrelated_split_witness
    @test Suslin.verify_sl3_local_murthy_input_context(split_direct_context)
    @test unrelated_split_source.target != ordinary.target
    wrong_split_witness_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        ordinary.target,
        ordinary.variable;
        split_witness = unrelated_split_witness,
    ))
    @test wrong_split_witness_err isa ArgumentError
    @test occursin(
        "split witness does not reconstruct the target",
        sprint(showerror, wrong_split_witness_err),
    )

    local_q0_unit = by_id["mg-local-q0-unit-at-u"]
    local_witness = first(local_q0_unit.witnesses)
    local_context = Suslin.sl3_local_murthy_input_context(
        local_q0_unit.entries.p,
        local_q0_unit.entries.q,
        local_q0_unit.entries.r,
        local_q0_unit.entries.s,
        local_q0_unit.variable;
        witness = local_witness,
    )
    @test local_context.target == local_q0_unit.target
    @test local_context.p0 == local_witness.p0
    @test local_context.q0 == local_witness.q0
    @test local_context.global_units.q0 == false
    @test local_context.local_units.q0 == true
    @test local_context.local_unit_witnesses.q0 == local_witness.local_unit_witness
    @test Suslin.verify_sl3_local_murthy_input_context(local_context)

    missing_local_witness_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        local_q0_unit.variable,
    ))
    @test missing_local_witness_err isa ArgumentError
    @test occursin("local-unit witness", sprint(showerror, missing_local_witness_err))

    wrong_variable_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        first(gens(local_q0_unit.ring.object));
        witness = local_witness,
    ))
    @test wrong_variable_err isa ArgumentError
    @test occursin("p must be monic", sprint(showerror, wrong_variable_err))

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    nonmonic_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        2 * X + one(R),
        X,
        R(2),
        one(R),
        X,
    ))
    @test nonmonic_err isa ArgumentError
    @test occursin("p must be monic", sprint(showerror, nonmonic_err))

    corrupted_q0_context = Suslin.SL3LocalMurthyInputContext(
        values(merge(
            _context_as_namedtuple(local_context),
            (; q0 = local_context.q0 + one(local_context.R)),
        ))...,
    )
    @test !Suslin.verify_sl3_local_murthy_input_context(corrupted_q0_context)

    corrupted_local_witness = merge(
        local_witness,
        (;
            local_unit_witness = merge(
                local_witness.local_unit_witness,
                (; residue_difference_coefficients = (zero(local_context.R),)),
            ),
        ),
    )
    @test_throws ArgumentError Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        local_q0_unit.variable;
        witness = corrupted_local_witness,
    )

    local_nonunit = by_id["mg-local-q0-nonunit-bezout-at-u"]
    bezout_witness = first(local_nonunit.witnesses)
    bezout_context = Suslin.sl3_local_murthy_input_context(
        local_nonunit.target,
        local_nonunit.variable;
        witness = bezout_witness,
    )
    @test bezout_context.global_units.branch_unit == false
    @test bezout_context.local_units.branch_unit == true
    @test bezout_context.bezout_witness == bezout_witness
    @test Suslin.verify_sl3_local_murthy_input_context(bezout_context)

    RU = local_nonunit.ring.object
    u = first(gens(RU))
    local_resultant_witness = _local_unit_witness_like(
        bezout_witness.branch_unit_witness,
        one(RU) + u,
        (one(RU),),
    )
    normalized_resultant_witnesses = Suslin._sl3_local_murthy_normalize_witness_data(
        (; resultant_unit_witness = local_resultant_witness),
        (;),
        nothing,
        nothing,
    )
    @test normalized_resultant_witnesses.local_unit_witnesses.resultant ==
        local_resultant_witness

    missing_resultant_unit_witness = (;
        p0 = bezout_witness.p0,
        q0 = bezout_witness.q0,
        p_prime = u * bezout_witness.p_prime,
        q_prime = u * bezout_witness.q_prime,
        resultant = u,
        p_prime_degree = 0,
        q_prime_degree = 1,
        branch_unit = 2 * u,
        case1_entries = (;
            p = local_nonunit.entries.p + u * bezout_witness.q_prime,
            q = local_nonunit.entries.q + u * bezout_witness.p_prime,
            r = u * bezout_witness.q_prime,
            s = u * bezout_witness.p_prime,
        ),
    )
    missing_resultant_err = _captured_error(
        () -> Suslin._sl3_local_murthy_validate_required_local_evidence(
            RU,
            2,
            1,
            u,
            (; q0 = false, resultant = false, branch_unit = false),
            (; q0 = false, resultant = false, branch_unit = false),
            missing_resultant_unit_witness,
        ),
    )
    @test missing_resultant_err isa ArgumentError
    @test occursin("Bezout resultant evidence", sprint(showerror, missing_resultant_err))

    missing_branch_unit_witness = (;
        p0 = bezout_witness.p0,
        q0 = bezout_witness.q0,
        p_prime = bezout_witness.p_prime,
        q_prime = bezout_witness.q_prime,
        resultant = bezout_witness.resultant,
        p_prime_degree = bezout_witness.p_prime_degree,
        q_prime_degree = bezout_witness.q_prime_degree,
        branch_unit = bezout_witness.branch_unit,
        case1_entries = bezout_witness.case1_entries,
    )
    missing_branch_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        local_nonunit.target,
        local_nonunit.variable;
        witness = missing_branch_unit_witness,
    ))
    @test missing_branch_err isa ArgumentError
    @test occursin("q0 or branch_unit", sprint(showerror, missing_branch_err))

    shifted_bezout_witness = merge(
        missing_branch_unit_witness,
        (;
            p_prime = bezout_witness.p_prime + local_nonunit.entries.q,
            q_prime = bezout_witness.q_prime + local_nonunit.entries.p,
            p_prime_degree = 1,
            q_prime_degree = 2,
            branch_unit = bezout_witness.branch_unit,
            case1_entries = (;
                p = local_nonunit.entries.p + bezout_witness.q_prime + local_nonunit.entries.p,
                q = local_nonunit.entries.q + bezout_witness.p_prime + local_nonunit.entries.q,
                r = bezout_witness.q_prime + local_nonunit.entries.p,
                s = bezout_witness.p_prime + local_nonunit.entries.q,
            ),
        ),
    )
    degree_guard_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        local_nonunit.target,
        local_nonunit.variable;
        witness = shifted_bezout_witness,
    ))
    @test degree_guard_err isa ArgumentError
    @test occursin("degree guards", sprint(showerror, degree_guard_err))

    Runit, (T,) = Oscar.polynomial_ring(QQ, ["T"])
    @test Suslin._sl3_local_murthy_bezout_data(
        Runit,
        T,
        one(Runit),
        zero(Runit),
        one(Runit),
        zero(Runit),
        0,
        -1,
        nothing,
    ) === nothing

    corrupted_bezout_witness = merge(
        bezout_witness,
        (; p_prime = bezout_witness.p_prime + one(local_nonunit.ring.object)),
    )
    @test_throws ArgumentError Suslin.sl3_local_murthy_input_context(
        local_nonunit.target,
        local_nonunit.variable;
        witness = corrupted_bezout_witness,
    )
end
