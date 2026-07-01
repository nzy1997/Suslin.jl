using Test
using Oscar
using Suslin

const ECP_COLUMN_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl")
include(ECP_COLUMN_CATALOG_PATH)

function _case_by_id(id::AbstractString)
    return ECPColumnFixtureCatalog.cases_by_id()[id]
end

function _column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _gf2_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    v = _column(entry)
    G = y * v[2] + v[3]
    return (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((; probe_id = :gf2_fixture_probe, G, lifted_tail_coefficients = (y, one(R)), tilde_G = G),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
end

function _qq_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    return (;
        source = :supplied_link_witness,
        residue_probes = (
            (; id = :qq_y_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),
            (; id = :qq_x_probe, kind = :deterministic_fixture, maximal_ideal_generators = (x,)),
        ),
        tail_reductions = (
            (; probe_id = :qq_y_probe, G = y, lifted_tail_coefficients = (zero(R), one(R)), tilde_G = y),
            (; probe_id = :qq_x_probe, G = x, lifted_tail_coefficients = (one(R), zero(R)), tilde_G = x),
        ),
        resultants = (y^2, y + one(R)),
        bezout_coefficients = ((; f = zero(R), h = y), (; f = one(R), h = -x)),
        coverage_multipliers = (one(R), one(R) - y),
        path_points = (zero(R), y^2 * x, x),
    )
end

function _factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _apply_factors(factors, column, R)
    return _factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _normality_witness(lower_certificate, entry)
    R = lower_certificate.ring
    n = length(lower_certificate.original_column)
    lower_product = _factor_product(lower_certificate.factors, R, n)
    return (;
        source = :supplied_normality_witness,
        conjugator = inv(lower_product),
        sl2_indices = (n, 1),
        sl2_entry = entry,
    )
end

function _replace_record_field(record, field::Symbol, value)
    fields = fieldnames(typeof(record))
    values = [getfield(record, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown record field: $(field)")
    values[idx] = value
    return typeof(record)(values...)
end

function _replace_rewrite_field(rewrite, field::Symbol, value)
    haskey(rewrite, field) || error("unknown rewrite field: $(field)")
    return merge(rewrite, NamedTuple{(field,)}((value,)))
end

function _replace_conjugated_certificate_field(cert, field::Symbol, value)
    fields = fieldnames(typeof(cert))
    values = [getfield(cert, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown conjugated certificate field: $(field)")
    values[idx] = value
    return typeof(cert)(values...)
end

function _tamper_nested_first_factor(cert, R, n::Int)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(R, n)
    return _replace_conjugated_certificate_field(cert, :factors, factors)
end

function _fixture_certificate(id::AbstractString, witness_builder, normality_entry_builder)
    entry = _case_by_id(id)
    R = entry.ring.object
    column = _column(entry)
    link = Suslin.ecp_link_step_certificate(
        column,
        R;
        variable_order = entry.ring.generators,
        selected_variable = entry.ring.generators[1],
        supplied_link_witness = witness_builder(entry),
    )
    lower = Suslin.ecp_column_reduction_certificate(collect(link.lower_variable_column), R)
    witness = _normality_witness(lower, normality_entry_builder(R, entry.ring.generators))
    certificate = Suslin.ecp_induction_normality_certificate(
        column,
        R;
        link_step = link,
        lower_reduction = lower,
        normality_witness = witness,
    )
    return (; entry, R, column, link, lower, witness, certificate)
end

@testset "ECP induction and normality replay" begin
    gf2 = _fixture_certificate(
        "ecp-variable-change-monic-gf2",
        _gf2_link_witness,
        (R, generators) -> generators[2] + one(R),
    )
    @test Suslin.verify_ecp_link_step_certificate(gf2.link)
    @test Suslin.verify_ecp_column_reduction(gf2.lower)
    @test gf2.certificate.verification.overall_ok
    @test gf2.certificate.normality_rewrite.sl2_block != identity_matrix(gf2.R, 2)
    @test gf2.certificate.normality_rewrite.fixed_lower_column_ok
    @test _apply_factors(gf2.certificate.final_factors, gf2.column, gf2.R) ==
          Suslin._target_reduced_column(gf2.R, length(gf2.column))
    @test Suslin.verify_ecp_induction_normality_certificate(gf2.certificate)

    qq = _fixture_certificate(
        "ecp-monic-first-entry-qq",
        _qq_link_witness,
        (R, generators) -> generators[1] + generators[2],
    )
    @test Suslin.verify_ecp_link_step_certificate(qq.link)
    @test Suslin.verify_ecp_column_reduction(qq.lower)
    @test qq.certificate.verification.overall_ok
    @test qq.certificate.normality_rewrite.sl2_block != identity_matrix(qq.R, 2)
    @test qq.certificate.normality_rewrite.rewrite_product == _factor_product(
        qq.certificate.normality_rewrite.rewrite_factors,
        qq.R,
        length(qq.column),
    )
    @test qq.certificate.normality_certificate isa Suslin.ConjugatedElementaryNormalityCertificate
    @test qq.certificate.normality_rewrite.normality_certificate == qq.certificate.normality_certificate
    @test Suslin.verify_conjugate_elementary_certificate(qq.certificate.normality_certificate)
    @test qq.certificate.normality_rewrite.rewrite_factors == qq.certificate.normality_certificate.factors
    @test qq.certificate.normality_rewrite.rewrite_product == qq.certificate.normality_certificate.product
    @test _apply_factors(qq.certificate.final_factors, qq.column, qq.R) ==
          Suslin._target_reduced_column(qq.R, length(qq.column))
    @test Suslin.verify_ecp_induction_normality_certificate(qq.certificate)

    explicit_sequence = Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = copy(qq.lower.factors),
        normality_witness = qq.witness,
    )
    @test Suslin.verify_ecp_induction_normality_certificate(explicit_sequence)

    supplied_nested = Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
        normality_witness = qq.witness,
        normality_certificate = qq.certificate.normality_certificate,
    )
    @test Suslin.verify_ecp_induction_normality_certificate(supplied_nested)
    @test supplied_nested.normality_certificate == qq.certificate.normality_certificate
    @test supplied_nested.final_factors == qq.certificate.final_factors

    mismatched_nested = Suslin.realize_conjugate_elementary_certificate(
        qq.witness.conjugator,
        length(qq.column),
        1,
        qq.witness.sl2_entry + one(qq.R),
    )
    @test_throws ArgumentError Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
        normality_witness = qq.witness,
        normality_certificate = mismatched_nested,
    )

    implicit_lower = Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        normality_witness = qq.witness,
    )
    @test Suslin.verify_ecp_induction_normality_certificate(implicit_lower)

    reordered_witness = (;
        sl2_entry = qq.witness.sl2_entry,
        sl2_indices = qq.witness.sl2_indices,
        conjugator = qq.witness.conjugator,
        source = qq.witness.source,
    )
    reordered_witness_certificate = Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
        normality_witness = reordered_witness,
    )
    @test Suslin.verify_ecp_induction_normality_certificate(reordered_witness_certificate)

    tampered_lifted = copy(qq.certificate.lifted_lower_variable_factors)
    tampered_lifted[1] = identity_matrix(qq.R, length(qq.column))
    tampered_lifted_certificate = _replace_record_field(
        qq.certificate,
        :lifted_lower_variable_factors,
        tampered_lifted,
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_lifted_certificate)

    tampered_rewrite_factors = copy(qq.certificate.normality_rewrite.rewrite_factors)
    tampered_rewrite_factors[1] = identity_matrix(qq.R, length(qq.column))
    tampered_rewrite = _replace_rewrite_field(
        qq.certificate.normality_rewrite,
        :rewrite_factors,
        tampered_rewrite_factors,
    )
    tampered_rewrite_certificate = _replace_record_field(
        qq.certificate,
        :normality_rewrite,
        tampered_rewrite,
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_rewrite_certificate)

    tampered_nested_factor = _replace_record_field(
        qq.certificate,
        :normality_certificate,
        _tamper_nested_first_factor(qq.certificate.normality_certificate, qq.R, length(qq.column)),
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_nested_factor)

    tampered_rewrite_nested = _replace_rewrite_field(
        qq.certificate.normality_rewrite,
        :normality_certificate,
        _tamper_nested_first_factor(qq.certificate.normality_certificate, qq.R, length(qq.column)),
    )
    tampered_rewrite_nested_certificate = _replace_record_field(
        qq.certificate,
        :normality_rewrite,
        tampered_rewrite_nested,
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_rewrite_nested_certificate)

    bad_target = copy(qq.certificate.normality_certificate.conjugation_target)
    bad_target[1, 1] += one(qq.R)
    tampered_nested_target = _replace_record_field(
        qq.certificate,
        :normality_certificate,
        _replace_conjugated_certificate_field(
            qq.certificate.normality_certificate,
            :conjugation_target,
            bad_target,
        ),
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_nested_target)

    tampered_lower_sequence_certificate = _replace_record_field(
        qq.certificate,
        :lower_reduction_certificate,
        nothing,
    )
    tampered_lower_sequence_certificate = _replace_record_field(
        tampered_lower_sequence_certificate,
        :lower_variable_factors,
        [:not_a_matrix],
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_lower_sequence_certificate)

    tampered_shape_witness = (;
        source = :supplied_normality_witness,
        conjugator = identity_matrix(qq.R, 2),
        sl2_indices = (2, 1),
        sl2_entry = one(qq.R),
    )
    @test_throws ArgumentError Suslin._ecp_induction_normality_rewrite(
        tampered_shape_witness,
        [one(qq.R), zero(qq.R)],
        [],
        qq.R,
    )

    @test !Suslin.verify_ecp_induction_normality_certificate((;))

    identity_witness = merge(qq.witness, (; sl2_entry = zero(qq.R)))
    @test_throws ArgumentError Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
        normality_witness = identity_witness,
    )

    automatic_normality = Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
    )
    @test automatic_normality.normality_witness.source == :constructed_normality_witness
    @test Suslin.verify_ecp_induction_normality_certificate(automatic_normality)

    wrong_conjugator_witness = merge(qq.witness, (; conjugator = identity_matrix(qq.R, length(qq.column))))
    @test_throws ArgumentError Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
        normality_witness = wrong_conjugator_witness,
    )

    @test_throws ArgumentError Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = (; factors = qq.lower.factors),
        normality_witness = qq.witness,
    )

    non_elementary_lower_reduction = [_factor_product(qq.lower.factors, qq.R, length(qq.column))]
    @test_throws ArgumentError Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = non_elementary_lower_reduction,
        normality_witness = qq.witness,
    )
end
