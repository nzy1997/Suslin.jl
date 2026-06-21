using Test
using Suslin
using Oscar

const SL3_MURTHY_GUPTA_ACCEPTANCE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _mg_acceptance_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _mg_acceptance_degree(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _mg_acceptance_special_form_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _mg_acceptance_assert_elementary_sequence(target, factors)
    R = base_ring(target)
    @test !isempty(factors)
    for factor in factors
        @test nrows(factor) == 3
        @test ncols(factor) == 3
        @test Suslin._same_base_ring(base_ring(factor), R)
        nonzero_offdiagonal = 0
        for i in 1:3, j in 1:3
            if i == j
                @test factor[i, j] == one(R)
            elseif factor[i, j] != zero(R)
                nonzero_offdiagonal += 1
            end
        end
        @test nonzero_offdiagonal <= 1
    end
    @test _mg_acceptance_product(factors, R) == target
    @test Suslin.verify_factorization(target, factors)
end

function _mg_acceptance_pre_murthy_open_or_unit_pivot(entry)
    R = entry.ring.object
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    if s == one(R) && p == one(R) + q * r
        return :open_s_one
    elseif p == one(R) && s == one(R) + q * r
        return :open_p_one
    elseif is_unit(s)
        return :s_unit
    elseif is_unit(p)
        return :p_unit
    end
    throw(ArgumentError("pre-#61 open-slice/unit-pivot solver does not support this target"))
end

function _mg_acceptance_assert_supported_case(entry; kwargs...)
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    X = entry.variable
    target = entry.target
    R = entry.ring.object

    @test target == _mg_acceptance_special_form_target(R, p, q, r, s)
    @test det(target) == one(R)
    @test Suslin._is_monic_in_variable(p, findfirst(isequal(X), collect(gens(R))), R)
    @test !is_unit(p)
    @test !is_unit(s)

    certificate_from_matrix = Suslin.realize_sl3_local_certificate(target, X; kwargs...)
    certificate_from_entries = Suslin.realize_sl3_local_certificate(p, q, r, s, X; kwargs...)
    @test certificate_from_matrix.target == target
    @test certificate_from_entries.target == target
    @test certificate_from_matrix.branch == certificate_from_entries.branch
    @test Suslin.verify_sl3_local_realization(certificate_from_matrix)
    @test Suslin.verify_sl3_local_realization(certificate_from_entries)
    _mg_acceptance_assert_elementary_sequence(target, certificate_from_matrix.factors)

    factors_from_matrix = Suslin.realize_sl3_local(target, X; kwargs...)
    factors_from_entries = Suslin.realize_sl3_local(p, q, r, s, X; kwargs...)
    @test factors_from_matrix == certificate_from_matrix.factors
    @test factors_from_entries == certificate_from_entries.factors
    _mg_acceptance_assert_elementary_sequence(target, factors_from_matrix)
    return certificate_from_matrix
end

function _mg_acceptance_assert_q0_unit_certificate(certificate)
    @test certificate.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(certificate)

    if certificate.witness.normalization !== nothing
        normalization = certificate.witness.normalization
        @test Suslin.verify_sl3_local_q_degree_normalization(normalization)
        @test normalization.target == certificate.target
        @test normalization.selected_variable == certificate.selected_variable
        @test certificate.witness.normalized_certificate !== nothing
        @test certificate.witness.reduction === nothing
        _mg_acceptance_assert_q0_unit_certificate(certificate.witness.normalized_certificate)
        return nothing
    end

    reduction = certificate.witness.reduction
    @test reduction !== nothing
    @test Suslin.verify_sl3_local_murthy_q_unit_reduction(reduction)
    @test reduction.target == certificate.target
    @test reduction.selected_variable == certificate.selected_variable
    @test reduction.q0 * reduction.q0_inverse == one(base_ring(certificate.target))
    @test _mg_acceptance_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_p
    @test reduction.split_certificate.branch == :murthy_split_lemma
    @test Suslin.verify_sl3_local_realization(reduction.split_certificate)
    @test Suslin.verify_sl3_local_split_lemma_replay(reduction.split_certificate.witness.split)
    @test reduction.split_certificate.witness.split.split_id == :murthy_q0_unit_split
    return nothing
end

function _mg_acceptance_assert_resultant_certificate(certificate; expected_source::Symbol)
    @test certificate.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(certificate)
    reduction = certificate.witness.reduction
    @test reduction.witness_source == expected_source
    @test Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(reduction)
    @test reduction.resultant == one(base_ring(certificate.target))
    @test reduction.p_prime * certificate.target[1, 1] - reduction.q_prime * certificate.target[1, 2] ==
        one(base_ring(certificate.target))
    @test reduction.target == reduction.left_factor * reduction.bezout_target
    @test reduction.bezout_target == reduction.first_elementary_factor * reduction.child_link_target
    @test reduction.branch_unit * reduction.branch_unit_inverse == one(base_ring(certificate.target))
    @test _mg_acceptance_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_q
    @test _mg_acceptance_degree(reduction.q_prime, reduction.selected_variable) < reduction.degree_p
    _mg_acceptance_assert_q0_unit_certificate(reduction.child_certificate)
    return nothing
end

@testset "Issue 61 Murthy-Gupta local SL3 acceptance" begin
    include(SL3_MURTHY_GUPTA_ACCEPTANCE_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    normalization_entry = by_id["mg-q-degree-normalization"]
    q0_unit_entry = by_id["mg-q0-unit-recursion"]
    supplied_entry = by_id["mg-q0-nonunit-normalized-bezout-resultant"]
    extracted_entry = by_id["mg-q0-nonunit-extracted-bezout-resultant"]

    acceptance_entries = (normalization_entry, q0_unit_entry, supplied_entry, extracted_entry)
    @test length(acceptance_entries) >= 3
    @test_throws ArgumentError _mg_acceptance_pre_murthy_open_or_unit_pivot(normalization_entry)
    @test_throws ArgumentError _mg_acceptance_pre_murthy_open_or_unit_pivot(supplied_entry)

    normalization_certificate = _mg_acceptance_assert_supported_case(normalization_entry)
    @test normalization_certificate.branch == :murthy_q0_unit
    @test normalization_certificate.witness.normalization !== nothing
    _mg_acceptance_assert_q0_unit_certificate(normalization_certificate)

    q0_unit_certificate = _mg_acceptance_assert_supported_case(q0_unit_entry)
    @test q0_unit_certificate.witness.normalization === nothing
    _mg_acceptance_assert_q0_unit_certificate(q0_unit_certificate)

    supplied_certificate = _mg_acceptance_assert_supported_case(
        supplied_entry;
        murthy_q0_nonunit_witness = first(supplied_entry.witnesses),
    )
    _mg_acceptance_assert_resultant_certificate(
        supplied_certificate;
        expected_source = :supplied_bezout_witness,
    )

    extracted_certificate = _mg_acceptance_assert_supported_case(extracted_entry)
    _mg_acceptance_assert_resultant_certificate(
        extracted_certificate;
        expected_source = :extracted_bezout_witness,
    )
end

@testset "Issue 61 staged local SL3 unsupported boundary" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    p = 2 * X + one(R)
    q = X
    r = R(2)
    s = one(R)
    target = _mg_acceptance_special_form_target(R, p, q, r, s)
    @test det(target) == one(R)
    @test !Suslin._is_monic_in_variable(p, 1, R)

    err = try
        Suslin.realize_sl3_local(target, X)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("staged local SL_3 solver failure", sprint(showerror, err))
    @test occursin("p must be monic in X", sprint(showerror, err))
end
