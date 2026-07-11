using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl"))

function _bridge_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _bridge_apply_factors(factors, column, R)
    return _bridge_factor_product(factors, R, length(column)) *
           matrix(R, length(column), 1, collect(column))
end

function _bridge_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _bridge_tamper(certificate; kwargs...)
    fields = propertynames(certificate)
    values = map(field -> get(kwargs, field, getproperty(certificate, field)), fields)
    return Suslin.LaurentToPolynomialECPBridgeCertificate(values...)
end

function _bridge_tamper_conversion(certificate; kwargs...)
    fields = propertynames(certificate)
    values = map(field -> get(kwargs, field, getproperty(certificate, field)), fields)
    return Suslin.LaurentToPolynomialColumnCertificate(values...)
end

function _bridge_tamper_noether(certificate; kwargs...)
    fields = propertynames(certificate)
    values = map(field -> get(kwargs, field, getproperty(certificate, field)), fields)
    return Suslin.LaurentNoetherCertificate(values...)
end

function _bridge_replace_substitution_value(substitution_map, idx::Int, value)
    replacement = collect(substitution_map)
    replacement[idx] = merge(replacement[idx], (; value))
    return tuple(replacement...)
end

function _bridge_tamper_noether_inverse(conversion_certificate)
    R = conversion_certificate.ring
    noether = conversion_certificate.noether_certificate
    idx = noether.other_generator_index
    inverse_substitution = _bridge_replace_substitution_value(
        noether.inverse_substitution,
        idx,
        noether.inverse_substitution[idx].value + one(R),
    )
    bad_noether = _bridge_tamper_noether(noether; inverse_substitution)
    return _bridge_tamper_conversion(
        conversion_certificate;
        noether_certificate = bad_noether,
    )
end

function _bridge_modify_first_elementary_factor(factors, R)
    modified = copy(factors)
    idx = findfirst(
        factor -> Suslin._canonical_elementary_factor_record(factor).kind == :elementary,
        modified,
    )
    idx === nothing && error("expected an elementary factor to modify")
    record = Suslin._canonical_elementary_factor_record(modified[idx])
    modified[idx] = elementary_matrix(
        record.n,
        record.row,
        record.col,
        record.coefficient + one(R),
        R,
    )
    return modified
end

function _bridge_forged_successful_child(child)
    return Suslin.ECPColumnReductionCertificate(
        child.original_column,
        child.ring,
        child.stages,
        Any[],
        child.final_column,
        child.verification,
    )
end

function _bridge_drop_last(factors)
    @test !isempty(factors)
    return factors[1:(end - 1)]
end

function _bridge_conversion_certificate(column, selected_entry_index::Int, selected_generator)
    noether = Suslin._laurent_noether_certificate(column, selected_entry_index, selected_generator)
    return Suslin._laurent_to_polynomial_certificate(
        column,
        noether,
        selected_entry_index,
        selected_generator,
    )
end

function _bridge_certificate(column, selected_entry_index::Int, selected_generator)
    conversion_certificate =
        _bridge_conversion_certificate(column, selected_entry_index, selected_generator)
    return Suslin._laurent_to_polynomial_ecp_bridge_certificate(conversion_certificate)
end

function _bridge_assert_common_certificate(certificate)
    conversion = certificate.conversion_certificate
    R = conversion.ring
    P = conversion.polynomial_ring
    n = length(conversion.original_column)

    @test certificate.validation_status == :ok
    @test Suslin._validate_laurent_to_polynomial_ecp_bridge_certificate(certificate) == :ok
    @test Suslin._validate_laurent_to_polynomial_certificate(conversion) == :ok
    @test Suslin.verify_ecp_column_reduction(certificate.ordinary_child_certificate)
    @test certificate.ordinary_child_certificate.original_column ==
          collect(conversion.polynomial_column)
    @test certificate.ordinary_child_certificate.ring === P
    @test certificate.ordinary_factors == certificate.ordinary_child_certificate.factors

    expected_lifted_factors = [
        Suslin._laurent_to_polynomial_lift_factor(factor, conversion.factor_lift_metadata)
        for factor in certificate.ordinary_child_certificate.factors
    ]
    expected_inverse_lifted_factors = Suslin._ecp_substitute_factor_sequence(
        expected_lifted_factors,
        conversion.noether_certificate.inverse_substitution,
        R,
    )
    expected_inverse_conversion_factors = Suslin._ecp_substitute_factor_sequence(
        conversion.conversion_factors,
        conversion.noether_certificate.inverse_substitution,
        R,
    )
    expected_complete_sequence = vcat(
        expected_inverse_lifted_factors,
        expected_inverse_conversion_factors,
    )

    @test certificate.raw_lifted_laurent_factors == expected_lifted_factors
    @test certificate.inverse_substituted_lifted_factors ==
          expected_inverse_lifted_factors
    @test certificate.laurent_conversion_factors == conversion.conversion_factors
    @test certificate.inverse_substituted_conversion_factors ==
          expected_inverse_conversion_factors
    @test certificate.complete_factor_sequence == expected_complete_sequence
    @test certificate.recomputed_product ==
          _bridge_factor_product(expected_complete_sequence, R, n)
    @test certificate.target_basis_column == _bridge_target_column(R, n)
    @test _bridge_apply_factors(
        certificate.complete_factor_sequence,
        conversion.original_column,
        R,
    ) == certificate.target_basis_column
    @test certificate.replay_summary.overall_ok
    @test all(factor -> base_ring(factor) === P, certificate.ordinary_factors)
    @test all(factor -> base_ring(factor) === R, certificate.raw_lifted_laurent_factors)
    @test all(
        factor -> base_ring(factor) === R,
        certificate.inverse_substituted_lifted_factors,
    )
    @test all(factor -> base_ring(factor) === R, certificate.laurent_conversion_factors)
    @test all(
        factor -> base_ring(factor) === R,
        certificate.inverse_substituted_conversion_factors,
    )
end

@testset "Laurent-to-polynomial ECP bridge positive controls" begin
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    polynomial_column = [x, y, one(R) + x]
    polynomial_bridge = _bridge_certificate(polynomial_column, 1, x)
    _bridge_assert_common_certificate(polynomial_bridge)
    @test polynomial_bridge.conversion_certificate.original_column ==
          tuple(polynomial_column...)
    @test polynomial_bridge.conversion_certificate.selected_entry_index == 1
    @test polynomial_bridge.conversion_certificate.selected_generator == x
    @test !isempty(polynomial_bridge.ordinary_factors)
    @test !isempty(polynomial_bridge.raw_lifted_laurent_factors)
    @test !isempty(polynomial_bridge.inverse_substituted_lifted_factors)
    @test isempty(polynomial_bridge.laurent_conversion_factors)
    @test isempty(polynomial_bridge.inverse_substituted_conversion_factors)

    U, (u, v) = Suslin.suslin_laurent_polynomial_ring(QQ, ["u", "v"])
    unit_column = [u^-1, v^-1 + u, one(U) + v]
    unit_bridge = _bridge_certificate(unit_column, 1, u)
    _bridge_assert_common_certificate(unit_bridge)
    @test unit_bridge.conversion_certificate.original_column == tuple(unit_column...)
    @test tuple(string.(unit_bridge.conversion_certificate.polynomial_column)...) ==
          ("0", "0", "1")
    @test !isempty(unit_bridge.laurent_conversion_factors)
    @test !isempty(unit_bridge.inverse_substituted_conversion_factors)
end

@testset "Laurent-to-polynomial ECP bridge negative controls" begin
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    polynomial_bridge = _bridge_certificate([x, y, one(R) + x], 1, x)
    @test length(polynomial_bridge.raw_lifted_laurent_factors) > 1
    @test Suslin._validate_laurent_to_polynomial_ecp_bridge_certificate(
        _bridge_tamper(
            polynomial_bridge;
            raw_lifted_laurent_factors =
                reverse(polynomial_bridge.raw_lifted_laurent_factors),
        ),
    ) != :ok

    tampered_conversion =
        _bridge_tamper_noether_inverse(polynomial_bridge.conversion_certificate)
    @test Suslin._validate_laurent_to_polynomial_ecp_bridge_certificate(
        _bridge_tamper(polynomial_bridge; conversion_certificate = tampered_conversion),
    ) != :ok

    changed_ordinary_factors = _bridge_modify_first_elementary_factor(
        polynomial_bridge.ordinary_factors,
        polynomial_bridge.conversion_certificate.polynomial_ring,
    )
    @test Suslin._validate_laurent_to_polynomial_ecp_bridge_certificate(
        _bridge_tamper(polynomial_bridge; ordinary_factors = changed_ordinary_factors),
    ) != :ok

    U, (u, v) = Suslin.suslin_laurent_polynomial_ring(QQ, ["u", "v"])
    unit_bridge = _bridge_certificate([u^-1, v^-1 + u, one(U) + v], 1, u)
    missing_conversion_factors = _bridge_drop_last(unit_bridge.laurent_conversion_factors)
    @test Suslin._validate_laurent_to_polynomial_ecp_bridge_certificate(
        _bridge_tamper(unit_bridge; laurent_conversion_factors = missing_conversion_factors),
    ) != :ok

    forged_child =
        _bridge_forged_successful_child(polynomial_bridge.ordinary_child_certificate)
    @test forged_child.verification.overall_ok
    @test !Suslin.verify_ecp_column_reduction(forged_child)
    @test Suslin._validate_laurent_to_polynomial_ecp_bridge_certificate(
        _bridge_tamper(polynomial_bridge; ordinary_child_certificate = forged_child),
    ) != :ok
end

@testset "Laurent-to-polynomial ECP bridge catalog probe" begin
    cases = LaurentFixtureCatalog.laurent_to_poly_route_cases_by_id()
    entry = cases["laurent-to-poly-general-ecp"]
    selected_generator = entry.ring.generators[1]
    conversion = _bridge_conversion_certificate(
        collect(entry.source_column),
        entry.selected_entry_index,
        selected_generator,
    )
    @test conversion.selected_entry_index == entry.selected_entry_index
    @test tuple(string.(conversion.polynomial_column)...) ==
          ("x^4*y + x", "x^2 + x + 1", "x^6*y^2 + x^4*y + 1")

    bridge = Suslin._laurent_to_polynomial_ecp_bridge_certificate(conversion)
    _bridge_assert_common_certificate(bridge)
    @test bridge.ordinary_child_certificate.stages[end].kind == :rank_one_normality_unit
    @test _bridge_apply_factors(
        bridge.complete_factor_sequence,
        conversion.original_column,
        conversion.ring,
    ) == bridge.target_basis_column
end
