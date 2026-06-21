using Test
using Suslin
using Oscar

const SL3_RESULTANT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _resultant_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _degree_in_variable(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _as_namedtuple(value)
    names = propertynames(value)
    return NamedTuple{names}(Tuple(getproperty(value, name) for name in names))
end

function _assert_resultant_certificate(cert; expected_source::Symbol)
    @test cert.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(cert)
    @test Suslin.verify_factorization(cert.target, cert.factors)
    @test _resultant_product(cert.factors, base_ring(cert.target)) == cert.target

    reduction = cert.witness.reduction
    @test reduction.witness_source == expected_source
    @test Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(reduction)
    @test reduction.child_certificate.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(reduction.child_certificate)
    @test reduction.bezout_target == reduction.first_elementary_factor * reduction.child_link_target
    @test reduction.target == reduction.left_factor * reduction.bezout_target
    @test reduction.branch_unit * reduction.branch_unit_inverse == one(base_ring(cert.target))
    @test _degree_in_variable(reduction.p_prime, reduction.selected_variable) < reduction.degree_q
    @test _degree_in_variable(reduction.q_prime, reduction.selected_variable) < reduction.degree_p
end

@testset "Murthy q(0)-nonunit Bezout/resultant branch for local SL3" begin
    include(SL3_RESULTANT_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    supplied_fixture = by_id["mg-q0-nonunit-normalized-bezout-resultant"]
    supplied_cert = Suslin.realize_sl3_local_certificate(
        supplied_fixture.entries.p,
        supplied_fixture.entries.q,
        supplied_fixture.entries.r,
        supplied_fixture.entries.s,
        supplied_fixture.variable;
        murthy_q0_nonunit_witness = first(supplied_fixture.witnesses),
    )
    _assert_resultant_certificate(supplied_cert; expected_source = :supplied_bezout_witness)

    extracted_fixture = by_id["mg-q0-nonunit-extracted-bezout-resultant"]
    extracted_cert = Suslin.realize_sl3_local_certificate(
        extracted_fixture.entries.p,
        extracted_fixture.entries.q,
        extracted_fixture.entries.r,
        extracted_fixture.entries.s,
        extracted_fixture.variable,
    )
    _assert_resultant_certificate(extracted_cert; expected_source = :extracted_bezout_witness)

    supplied_reduction = _as_namedtuple(supplied_cert.witness.reduction)
    corrupt_supplied_reduction = merge(
        supplied_reduction,
        (;
            p_prime = supplied_reduction.p_prime + one(base_ring(supplied_cert.target)),
        ),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_supplied_reduction)

    extracted_reduction = _as_namedtuple(extracted_cert.witness.reduction)
    tampered_bezout_reduction = merge(
        extracted_reduction,
        (;
            bezout_target = extracted_reduction.bezout_target + identity_matrix(base_ring(extracted_cert.target), 3),
        ),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(tampered_bezout_reduction)

    tampered_child_reduction = merge(
        extracted_reduction,
        (;
            child_certificate = merge(
                _as_namedtuple(extracted_reduction.child_certificate),
                (; target = extracted_reduction.child_certificate.target + identity_matrix(base_ring(extracted_cert.target), 3)),
            ),
        ),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(tampered_child_reduction)
end
