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
