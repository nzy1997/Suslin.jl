using Test
using Oscar
using Suslin

include(joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl"))

function _general_case_by_id(id::AbstractString)
    return ECPMainlineFixtureCatalog.cases_by_id()[id]
end

function _general_column(entry)
    return [getproperty(entry.column_entries, name) for name in entry.column_order]
end

function _general_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _apply_factors(factors, column, R)
    return _general_factor_product(factors, R, length(column)) *
           matrix(R, length(column), 1, collect(column))
end

function _replace_record_field(record, field::Symbol, value)
    fields = fieldnames(typeof(record))
    values = [getfield(record, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown record field: $(field)")
    values[idx] = value
    return typeof(record)(values...)
end

function _replace_namedtuple_field(record::NamedTuple, field::Symbol, value)
    haskey(record, field) || error("unknown named tuple field: $(field)")
    return merge(record, NamedTuple{(field,)}((value,)))
end

function _replace_nested_field(record, field::Symbol, value)
    record isa NamedTuple && return _replace_namedtuple_field(record, field, value)
    return _replace_record_field(record, field, value)
end

function _replace_conjugated_certificate_field(cert, field::Symbol, value)
    return _replace_record_field(cert, field, value)
end

function _tamper_nested_first_factor(cert, R, n::Int)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(R, n)
    return _replace_conjugated_certificate_field(cert, :factors, factors)
end

function _same_context_case()
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    column = [y + one(R), one(R), zero(R)]
    witness = (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :same_context_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((;
            probe_id = :same_context_probe,
            G = one(R),
            lifted_tail_coefficients = (one(R), zero(R)),
            tilde_G = one(R),
        ),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = zero(R), h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
    link = Suslin.ecp_link_step_certificate(
        column,
        R;
        variable_order = (x, y),
        selected_variable = x,
        supplied_link_witness = witness,
    )
    return (; R, column, link)
end

function _unsupported_lower_case()
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    column = [x + y, x^2 + one(R), x * y + one(R), zero(R)]
    witness = (;
        source = :supplied_link_witness,
        residue_probes = (
            (; id = :unsupported_y_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),
            (; id = :unsupported_x_probe, kind = :deterministic_fixture, maximal_ideal_generators = (x,)),
        ),
        tail_reductions = (
            (; probe_id = :unsupported_y_probe, G = x^2 + one(R), lifted_tail_coefficients = (one(R), zero(R), zero(R)), tilde_G = x^2 + one(R)),
            (; probe_id = :unsupported_x_probe, G = x * y + one(R), lifted_tail_coefficients = (zero(R), one(R), zero(R)), tilde_G = x * y + one(R)),
        ),
        resultants = (one(R), one(R)),
        bezout_coefficients = (
            (; f = zero(R), h = one(R)),
            (; f = zero(R), h = one(R)),
        ),
        coverage_multipliers = (one(R), zero(R)),
        path_points = (zero(R), x, x),
    )
    link_witness = Suslin.ecp_link_witness(
        column,
        R;
        variable_order = (x, y),
        selected_variable = x,
        supplied_link_witness = witness,
    )
    link = Suslin.ecp_link_step_certificate(column, R; link_witness = link_witness)
    return (; R, column, link)
end

@testset "General ECP induction/normality composition" begin
    entry = _general_case_by_id("ecp-mainline-sl3-route-qq")
    R = entry.ring.object
    column = _general_column(entry)
    witness = entry.support_evidence.link_witness
    link = Suslin.ecp_link_step_certificate(
        column,
        R;
        link_witness = witness,
        route_mode = :polynomial_sl3,
    )
    lower = Suslin.ecp_column_reduction_certificate(collect(link.lower_variable_column), R)

    cert = Suslin.ecp_induction_normality_certificate(
        column,
        R;
        link_step = link,
    )

    @test cert isa Suslin.ECPInductionNormalityCertificate
    @test cert.descent_measure.strict_descent
    @test cert.descent_measure.parent_profile != cert.descent_measure.lower_profile
    @test cert.lower_reduction_certificate isa Suslin.ECPColumnReductionCertificate
    @test Suslin.verify_ecp_column_reduction(cert.lower_reduction_certificate)
    @test cert.lower_reduction_certificate.original_column == collect(link.lower_variable_column)
    @test cert.normality_witness.source == :constructed_normality_witness
    @test cert.normality_certificate isa Suslin.ConjugatedElementaryNormalityCertificate
    @test Suslin.verify_conjugate_elementary_certificate(cert.normality_certificate)
    @test cert.normality_rewrite.normality_certificate == cert.normality_certificate
    @test Suslin.verify_ecp_induction_normality_certificate(cert)
    @test _apply_factors(cert.final_factors, column, R) ==
          Suslin._target_reduced_column(R, length(column))

    tampered_descent = _replace_record_field(
        cert,
        :descent_measure,
        _replace_nested_field(cert.descent_measure, :strict_descent, false),
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_descent)

    tampered_lower_certificate = _replace_record_field(cert, :lower_reduction_certificate, nothing)
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_lower_certificate)

    tampered_lifted_factors = copy(cert.lifted_lower_variable_factors)
    tampered_lifted_factors[1] = identity_matrix(R, length(column))
    @test !Suslin.verify_ecp_induction_normality_certificate(
        _replace_record_field(cert, :lifted_lower_variable_factors, tampered_lifted_factors),
    )

    tampered_normality_witness = _replace_record_field(
        cert,
        :normality_witness,
        merge(cert.normality_witness, (; sl2_entry = cert.normality_witness.sl2_entry + one(R))),
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_normality_witness)

    tampered_normality_certificate = _replace_record_field(
        cert,
        :normality_certificate,
        _tamper_nested_first_factor(cert.normality_certificate, R, length(column)),
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_normality_certificate)

    tampered_final_factors = copy(cert.final_factors)
    tampered_final_factors[1] = identity_matrix(R, length(column))
    @test !Suslin.verify_ecp_induction_normality_certificate(
        _replace_record_field(cert, :final_factors, tampered_final_factors),
    )

    same_context = _same_context_case()
    @test_throws Regex("same-context recursive lower-variable call") Suslin.ecp_induction_normality_certificate(
        same_context.column,
        same_context.R;
        link_step = same_context.link,
    )

    unsupported = _unsupported_lower_case()
    @test_throws Regex("missing lower-variable reduction") Suslin.ecp_induction_normality_certificate(
        unsupported.column,
        unsupported.R;
        link_step = unsupported.link,
    )

    valid_normality_witness = (;
        source = :supplied_normality_witness,
        conjugator = inv(_general_factor_product(lower.factors, R, length(column))),
        sl2_indices = (length(column), 1),
        sl2_entry = one(R),
    )
    mismatched_normality_certificate = Suslin.realize_conjugate_elementary_certificate(
        valid_normality_witness.conjugator,
        length(column),
        1,
        valid_normality_witness.sl2_entry + one(R),
    )
    @test_throws Regex("missing normality rewrite") Suslin.ecp_induction_normality_certificate(
        column,
        R;
        link_step = link,
        lower_reduction = lower,
        normality_witness = valid_normality_witness,
        normality_certificate = mismatched_normality_certificate,
    )
end
