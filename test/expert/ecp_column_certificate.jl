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

function _ecp_cert_monicity_certificate(column, R)
    result = Suslin._reduce_after_monicity_normalization_certificate(column, R)
    @test result !== nothing
    return Suslin._ecp_certificate_from_stage(column, R, result.stage)
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

    monic_column = _ecp_cert_column(cases["ecp-variable-change-monic-gf2"])
    monic_ring = cases["ecp-variable-change-monic-gf2"].ring.object
    monic_cert = Suslin.ecp_column_reduction_certificate(monic_column, monic_ring)
    @test any(stage -> stage.kind == :ecp_pipeline, monic_cert.stages)
    @test monic_cert.stages[end].route_metadata.route == :general_ecp_pipeline
    _assert_ecp_certificate_replays(monic_cert)
    monicity_cert = _ecp_cert_monicity_certificate(monic_column, monic_ring)
    @test any(stage -> stage.kind == :monicity_normalization, monicity_cert.stages)

    S, (t,) = Oscar.polynomial_ring(GF(2), ["t"])
    embedded_column = [zero(S), t^2, t^3, t^2 + t + one(S)]
    embedded_cert = Suslin.ecp_column_reduction_certificate(embedded_column, S)
    @test any(stage -> stage.kind == :embedded_three_block, embedded_cert.stages)
    _assert_ecp_certificate_replays(embedded_cert)

    legacy_embedded_column = [one(S), t, t^2, t + one(S)]
    legacy_subfactors = Suslin._reduce_exact_small_column(legacy_embedded_column[1:3], S)
    legacy_embedded_factors = Suslin._embedded_three_block_reduction(
        legacy_embedded_column,
        S,
        (1, 2, 3),
        legacy_subfactors,
    )
    @test _ecp_cert_apply_factors(legacy_embedded_factors, legacy_embedded_column, S) ==
          _ecp_cert_target_column(S, length(legacy_embedded_column))

    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    laurent_column = [
        x^-1 * y^-1 * (x + y^2),
        x^-1 * y^-1 * (x * y + x + one(R)),
        x^-1 * y^-1 * (x^2 + x * y + y + one(R)),
    ]
    laurent_cert = Suslin.ecp_column_reduction_certificate(laurent_column, R)
    @test any(stage -> stage.kind == :laurent_normalization, laurent_cert.stages)
    _assert_ecp_certificate_replays(laurent_cert)

    @test !Suslin.verify_ecp_column_reduction(_tamper_first_factor(unit_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_witness(witness_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_witness_with_extra_entry(witness_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_monicity_inverse(monicity_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_embedded_indices(embedded_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_laurent_shift(laurent_cert))
    @test !Suslin.verify_ecp_column_reduction((; original_column = unit_cert.original_column))

    supported_factors = Suslin._reduce_supported_unimodular_column(unit_cert.original_column, unit_cert.ring)
    @test _ecp_cert_apply_factors(supported_factors, unit_cert.original_column, unit_cert.ring) ==
          _ecp_cert_target_column(unit_cert.ring, length(unit_cert.original_column))

    witness_stage = only(stage for stage in witness_cert.stages if stage.kind == :witness_unit)
    witness_factors = Suslin._reduce_via_witness_unit(
        witness_cert.original_column,
        collect(witness_stage.witness),
        witness_stage.pivot_index,
        witness_cert.ring,
    )
    @test _ecp_cert_apply_factors(witness_factors, witness_cert.original_column, witness_cert.ring) ==
          _ecp_cert_target_column(witness_cert.ring, length(witness_cert.original_column))

    polynomial_factors = Suslin._reduce_polynomial_unimodular_column_exact(monic_cert.original_column, monic_cert.ring)
    @test _ecp_cert_apply_factors(polynomial_factors, monic_cert.original_column, monic_cert.ring) ==
          _ecp_cert_target_column(monic_cert.ring, length(monic_cert.original_column))

    exact_small_factors = Suslin._reduce_exact_small_column(monic_cert.original_column, monic_cert.ring)
    @test _ecp_cert_apply_factors(exact_small_factors, monic_cert.original_column, monic_cert.ring) ==
          _ecp_cert_target_column(monic_cert.ring, length(monic_cert.original_column))

    monicity_factors = Suslin._reduce_after_monicity_normalization(monic_cert.original_column, monic_cert.ring)
    @test _ecp_cert_apply_factors(monicity_factors, monic_cert.original_column, monic_cert.ring) ==
          _ecp_cert_target_column(monic_cert.ring, length(monic_cert.original_column))

    unknown_stage_cert = Suslin.ECPColumnReductionCertificate(
        unit_cert.original_column,
        unit_cert.ring,
        (unit_cert.stages[1], (; kind = :unknown_stage)),
        Any[],
        unit_cert.final_column,
        nothing,
    )
    @test !Suslin.verify_ecp_column_reduction(unknown_stage_cert)

    unsupported = _ecp_cert_column(cases["ecp-unsupported-unimodular-gf2"])
    @test_throws ArgumentError Suslin.ecp_column_reduction_certificate(unsupported, cases["ecp-unsupported-unimodular-gf2"].ring.object)

    non_unimodular = _ecp_cert_column(cases["ecp-non-unimodular-gf2"])
    @test_throws ArgumentError Suslin.ecp_column_reduction_certificate(non_unimodular, cases["ecp-non-unimodular-gf2"].ring.object)
end
