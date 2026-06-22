using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))

function _ecp_cert_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _ecp_cert_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _ecp_cert_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _ecp_cert_apply_factors(factors, column, R)
    return _ecp_cert_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _assert_ecp_certificate_replays(cert)
    R = cert.ring
    @test cert isa Suslin.ECPColumnReductionCertificate
    @test Suslin.verify_ecp_column_reduction(cert)
    @test _ecp_cert_apply_factors(cert.factors, cert.original_column, R) == _ecp_cert_target_column(R, length(cert.original_column))
    @test cert.final_column == _ecp_cert_target_column(R, length(cert.original_column))
    legacy_factors = Suslin.reduce_unimodular_column(cert.original_column, R)
    @test legacy_factors isa Vector
    @test _ecp_cert_apply_factors(legacy_factors, cert.original_column, R) == _ecp_cert_target_column(R, length(cert.original_column))
end

function _tamper_first_factor(cert)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(cert.ring, length(cert.original_column))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        cert.stages,
        factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_witness(cert)
    stages = collect(cert.stages)
    witness_idx = findfirst(stage -> stage.kind == :witness_unit, stages)
    witness_idx === nothing && error("certificate has no witness stage")
    stage = stages[witness_idx]
    witness = collect(stage.witness)
    witness[1] += one(cert.ring)
    stages[witness_idx] = merge(stage, (; witness = tuple(witness...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_witness_with_extra_entry(cert)
    stages = collect(cert.stages)
    witness_idx = findfirst(stage -> stage.kind == :witness_unit, stages)
    witness_idx === nothing && error("certificate has no witness stage")
    stage = stages[witness_idx]
    witness = tuple(stage.witness..., one(cert.ring))
    stages[witness_idx] = merge(stage, (; witness))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_monicity_inverse(cert)
    stages = collect(cert.stages)
    monic_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    monic_idx === nothing && error("certificate has no monicity stage")
    stage = stages[monic_idx]
    inverse_values = collect(stage.inverse_values)
    inverse_values[stage.variable_index] = inverse_values[stage.variable_index] + one(cert.ring)
    stages[monic_idx] = merge(stage, (; inverse_values = tuple(inverse_values...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_embedded_indices(cert)
    stages = collect(cert.stages)
    embedded_idx = findfirst(stage -> stage.kind == :embedded_three_block, stages)
    embedded_idx === nothing && error("certificate has no embedded block stage")
    stage = stages[embedded_idx]
    stages[embedded_idx] = merge(stage, (; indices = reverse(stage.indices)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_laurent_shift(cert)
    stages = collect(cert.stages)
    laurent_idx = findfirst(stage -> stage.kind == :laurent_normalization, stages)
    laurent_idx === nothing && error("certificate has no Laurent normalization stage")
    stage = stages[laurent_idx]
    stages[laurent_idx] = merge(stage, (; shift = one(cert.ring),))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

@testset "ECP column reduction certificates" begin
    cases = ECPColumnFixtureCatalog.cases_by_id()

    unit_cert = Suslin.ecp_column_reduction_certificate(_ecp_cert_column(cases["ecp-unit-entry-gf2"]), cases["ecp-unit-entry-gf2"].ring.object)
    @test any(stage -> stage.kind == :unit_entry, unit_cert.stages)
    _assert_ecp_certificate_replays(unit_cert)

    witness_cert = Suslin.ecp_column_reduction_certificate(_ecp_cert_column(cases["ecp-witness-unit-gf2"]), cases["ecp-witness-unit-gf2"].ring.object)
    @test any(stage -> stage.kind == :witness_unit, witness_cert.stages)
    _assert_ecp_certificate_replays(witness_cert)

    monic_cert = Suslin.ecp_column_reduction_certificate(_ecp_cert_column(cases["ecp-variable-change-monic-gf2"]), cases["ecp-variable-change-monic-gf2"].ring.object)
    @test any(stage -> stage.kind == :monicity_normalization, monic_cert.stages)
    _assert_ecp_certificate_replays(monic_cert)

    S, (t,) = Oscar.polynomial_ring(GF(2), ["t"])
    embedded_column = [zero(S), t^2, t^3, t^2 + t + one(S)]
    embedded_cert = Suslin.ecp_column_reduction_certificate(embedded_column, S)
    @test any(stage -> stage.kind == :embedded_three_block, embedded_cert.stages)
    _assert_ecp_certificate_replays(embedded_cert)

    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    laurent_column = [
        x^-1 + x^-2 * y^2,
        x^-1 * y + x^-1 + x^-2,
        one(R) + x^-1 * y + x^-2 * y + x^-2,
        x^-1 + x^-2 * y,
        x^-1 * y + x^-2 * y^2,
        x^-2 * y + x^-1 * y^2,
    ]
    laurent_cert = Suslin.ecp_column_reduction_certificate(laurent_column, R)
    @test any(stage -> stage.kind == :laurent_normalization, laurent_cert.stages)
    _assert_ecp_certificate_replays(laurent_cert)

    @test !Suslin.verify_ecp_column_reduction(_tamper_first_factor(unit_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_witness(witness_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_witness_with_extra_entry(witness_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_monicity_inverse(monic_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_embedded_indices(embedded_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_laurent_shift(laurent_cert))
    @test !Suslin.verify_ecp_column_reduction((; original_column = unit_cert.original_column))

    unsupported = _ecp_cert_column(cases["ecp-unsupported-unimodular-gf2"])
    @test_throws ArgumentError Suslin.ecp_column_reduction_certificate(unsupported, cases["ecp-unsupported-unimodular-gf2"].ring.object)

    non_unimodular = _ecp_cert_column(cases["ecp-non-unimodular-gf2"])
    @test_throws ArgumentError Suslin.ecp_column_reduction_certificate(non_unimodular, cases["ecp-non-unimodular-gf2"].ring.object)
end
