using Test
using Oscar
using Suslin

function _ecp_acceptance_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _ecp_acceptance_apply(factors, column, R)
    return _ecp_acceptance_product(factors, R, length(column)) *
        matrix(R, length(column), 1, collect(column))
end

function _ecp_acceptance_target(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _ecp_acceptance_assert_public(column, R)
    @test !any(is_unit, column)
    @test Suslin._reduce_supported_unimodular_column_certificate(column, R) === nothing
    factors = Suslin.reduce_unimodular_column(column, R)
    @test _ecp_acceptance_apply(factors, column, R) == _ecp_acceptance_target(R, length(column))
    return factors
end

function _ecp_acceptance_gf2_cases()
    R, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    base = [
        x + y^2,
        x * y + x + one(R),
        x^2 + x * y + y + one(R),
    ]
    return R, [
        ("canonical-full-route", base[[1, 2, 3]]),
        ("inverse-substitution-permuted", base[[2, 1, 3]]),
        ("permuted-third-middle", base[[1, 3, 2]]),
        ("permuted-cyclic", base[[3, 1, 2]]),
    ]
end

function _ecp_acceptance_length4_general_case()
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    column = [
        x * y,
        x * (one(R) - y),
        (one(R) - x) * y,
        (one(R) - x) * (one(R) - y),
    ]
    for (target, source, coeff) in (
        (2, 4, -x * y - y),
        (3, 2, -x * y - y),
        (4, 3, x^2),
        (2, 3, x * y + y),
        (2, 3, -y - one(R)),
    )
        column[target] += coeff * column[source]
    end
    return R, column
end

function _ecp_acceptance_good_link_witness(column, R)
    x, y = gens(R)
    G = y * column[2] + column[3]
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

function _ecp_acceptance_normality_witness(lower_certificate, entry)
    R = lower_certificate.ring
    n = length(lower_certificate.original_column)
    lower_product = _ecp_acceptance_product(lower_certificate.factors, R, n)
    return (;
        source = :supplied_normality_witness,
        conjugator = inv(lower_product),
        sl2_indices = (n, 1),
        sl2_entry = entry,
    )
end

function _ecp_acceptance_capture_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _ecp_acceptance_tampered_staged_certificate(
    staged;
    monicity = staged.monicity,
    lower_reduction = staged.lower_reduction,
    normality_witness = staged.normality_witness,
)
    return Suslin.ECPStagedColumnReductionCertificate(
        staged.original_column,
        staged.ring,
        monicity,
        staged.link_step,
        lower_reduction,
        normality_witness,
        staged.induction_normality,
        staged.factors,
        staged.final_column,
        staged.verification,
    )
end

@testset "public ECP unimodular-column pipeline" begin
    R, cases = _ecp_acceptance_gf2_cases()
    public_factors_by_name = Dict{String, Any}()
    for (name, column) in cases
        public_factors_by_name[name] = _ecp_acceptance_assert_public(column, R)
    end

    canonical = cases[1][2]
    staged = Suslin.ecp_staged_column_reduction_certificate(canonical, R)
    @test Suslin.verify_ecp_staged_column_reduction(staged)
    @test !Suslin.verify_ecp_staged_column_reduction((; malformed = true))
    @test staged.verification.route == (:validation, :monicity_forcing, :link_witness, :link_step, :induction_normality)
    @test staged.verification.link_witness_ok
    @test staged.verification.link_step_ok
    @test staged.verification.induction_normality_ok
    @test staged.induction_normality.normality_rewrite.sl2_block != identity_matrix(R, 2)
    @test staged.induction_normality.normality_certificate isa Suslin.ConjugatedElementaryNormalityCertificate
    @test staged.induction_normality.normality_rewrite.normality_certificate ==
          staged.induction_normality.normality_certificate
    @test Suslin.verify_conjugate_elementary_certificate(staged.induction_normality.normality_certificate)
    public_cert = Suslin.ecp_column_reduction_certificate(canonical, R)
    @test Suslin.verify_ecp_column_reduction(public_cert)
    @test public_factors_by_name["canonical-full-route"] == public_cert.factors
    @test public_cert.stages[end].kind == :ecp_pipeline
    @test public_cert.stages[end].route_metadata.route == :general_ecp_pipeline
    @test public_cert.factors != staged.factors
    canonical_diagnostic = Suslin.diagnose_unimodular_column_reduction(canonical, R)
    @test canonical_diagnostic.status == :supported
    @test canonical_diagnostic.attempted_stages[end] == :general_ecp_pipeline
    canonical_general_detail = only(filter(
        detail -> detail.stage == :general_ecp_pipeline,
        canonical_diagnostic.stage_details,
    ))
    @test canonical_general_detail.outcome == :supported
    @test canonical_general_detail.link_route_mode == public_cert.stages[end].route_metadata.link_route_mode

    general_R, general_column = _ecp_acceptance_length4_general_case()
    @test length(general_column) > 3
    @test Suslin.is_unimodular_column(general_column, general_R)
    @test !any(is_unit, general_column)
    @test Suslin._reduce_supported_unimodular_column_certificate(general_column, general_R) === nothing
    @test Suslin._reduce_via_supported_three_block_certificate(general_column, general_R) === nothing

    general_diagnostic = Suslin.diagnose_unimodular_column_reduction(general_column, general_R)
    @test general_diagnostic.status == :supported
    @test general_diagnostic.attempted_stages[end] == :general_ecp_pipeline
    general_detail = only(filter(
        detail -> detail.stage == :general_ecp_pipeline,
        general_diagnostic.stage_details,
    ))
    @test general_detail.outcome == :supported
    @test general_detail.link_route_mode == :direct_elementary
    @test general_detail.normalized_column_length == length(general_column)

    general_cert = Suslin.ecp_column_reduction_certificate(general_column, general_R)
    @test Suslin.verify_ecp_column_reduction(general_cert)
    @test general_cert.stages[end].kind == :ecp_pipeline
    @test general_cert.stages[end].route_metadata.route == :general_ecp_pipeline
    @test general_cert.stages[end].route_metadata.link_route_mode == :direct_elementary
    @test general_cert.stages[end].route_metadata.normalized_column_length == length(general_column)
    @test _ecp_acceptance_apply(general_cert.factors, general_column, general_R) ==
          _ecp_acceptance_target(general_R, length(general_column))
    @test Suslin.reduce_unimodular_column(general_column, general_R) == general_cert.factors

    tampered_factors = copy(general_cert.factors)
    tampered_factors[1] = identity_matrix(general_R, length(general_column))
    tampered_cert = Suslin.ECPColumnReductionCertificate(
        general_cert.original_column,
        general_cert.ring,
        general_cert.stages,
        tampered_factors,
        general_cert.final_column,
        general_cert.verification,
    )
    @test !Suslin.verify_ecp_column_reduction(tampered_cert)

    explicit_public_staged = Suslin._ecp_public_staged_reduction_certificate(
        canonical,
        R;
        supplied_link_witness = _ecp_acceptance_good_link_witness(canonical, R),
    )
    @test Suslin.verify_ecp_staged_column_reduction(explicit_public_staged)
    @test _ecp_acceptance_apply(explicit_public_staged.factors, canonical, R) ==
          _ecp_acceptance_target(R, length(canonical))

    tampered_monicity = _ecp_acceptance_tampered_staged_certificate(
        staged;
        monicity = merge(staged.monicity, (; source = :tampered_monicity)),
    )
    @test !Suslin.verify_ecp_staged_column_reduction(tampered_monicity)

    tampered_lower_reduction = _ecp_acceptance_tampered_staged_certificate(
        staged;
        lower_reduction = Any[],
    )
    @test !Suslin.verify_ecp_staged_column_reduction(tampered_lower_reduction)

    raw_lower_link = Suslin.ecp_link_step_certificate(
        canonical,
        R;
        variable_order = gens(R),
        selected_variable = first(gens(R)),
        supplied_link_witness = _ecp_acceptance_good_link_witness(canonical, R),
    )
    raw_lower = Suslin.ecp_column_reduction_certificate(collect(raw_lower_link.lower_variable_column), R)
    staged_with_raw_lower = Suslin.ecp_staged_column_reduction_certificate(
        canonical,
        R;
        supplied_link_witness = _ecp_acceptance_good_link_witness(canonical, R),
        lower_reduction = raw_lower.factors,
    )
    @test Suslin.verify_ecp_staged_column_reduction(staged_with_raw_lower)
    tampered_missing_raw_lower = _ecp_acceptance_tampered_staged_certificate(
        staged_with_raw_lower;
        lower_reduction = nothing,
    )
    @test !Suslin.verify_ecp_staged_column_reduction(tampered_missing_raw_lower)

    tampered_normality = _ecp_acceptance_tampered_staged_certificate(
        staged;
        normality_witness = merge(staged.normality_witness, (; source = :tampered_normality)),
    )
    @test !Suslin.verify_ecp_staged_column_reduction(tampered_normality)

    permuted = cases[2][2]
    permuted_cert = Suslin.ecp_column_reduction_certificate(permuted, R)
    pipeline_stage = only([stage for stage in permuted_cert.stages if stage.kind == :ecp_pipeline])
    @test pipeline_stage.route_metadata.route == :general_ecp_pipeline
    @test all(
        factor -> base_ring(factor) == R,
        pipeline_stage.normalization.inverse_substituted_coordinate_move_factors,
    )
    @test all(factor -> base_ring(factor) == R, pipeline_stage.inverse_substituted_induction_factors)
    @test _ecp_acceptance_apply(pipeline_stage.factors, permuted, R) ==
          _ecp_acceptance_target(R, length(permuted))

    x, y = gens(R)
    non_unimodular = [x, y, x * y]
    err = _ecp_acceptance_capture_error(() -> Suslin.reduce_unimodular_column(non_unimodular, R))
    @test err isa ArgumentError
    @test occursin("v must be a unimodular column", sprint(showerror, err))

    bad_link = merge(_ecp_acceptance_good_link_witness(canonical, R), (; resultants = (x,),))
    @test_throws ArgumentError Suslin.ecp_staged_column_reduction_certificate(
        canonical,
        R;
        supplied_link_witness = bad_link,
    )

    err = _ecp_acceptance_capture_error(
        () -> Suslin.ecp_staged_column_reduction_certificate(
            canonical,
            R;
            lower_reduction = Any[],
        ),
    )
    @test err isa ArgumentError
    @test (
        occursin("lower-variable reduction", sprint(showerror, err)) ||
        occursin("v(0)", sprint(showerror, err))
    )

    err = _ecp_acceptance_capture_error(
        () -> Suslin._ecp_public_staged_reduction_certificate(
            canonical,
            R;
            lower_reduction = Any[],
        ),
    )
    @test err isa ArgumentError
    @test (
        occursin("lower-variable reduction", sprint(showerror, err)) ||
        occursin("v(0)", sprint(showerror, err))
    )

    good_link = Suslin.ecp_link_step_certificate(
        canonical,
        R;
        variable_order = gens(R),
        selected_variable = x,
        supplied_link_witness = _ecp_acceptance_good_link_witness(canonical, R),
    )
    lower = Suslin.ecp_column_reduction_certificate(collect(good_link.lower_variable_column), R)
    bad_normality = merge(
        _ecp_acceptance_normality_witness(lower, y + one(R)),
        (; conjugator = identity_matrix(R, length(canonical))),
    )
    @test_throws ArgumentError Suslin.ecp_staged_column_reduction_certificate(
        canonical,
        R;
        supplied_link_witness = _ecp_acceptance_good_link_witness(canonical, R),
        normality_witness = bad_normality,
    )

    err = _ecp_acceptance_capture_error(
        () -> Suslin.ecp_staged_column_reduction_certificate(
            canonical,
            R;
            variable_order = (),
        ),
    )
    @test err isa ArgumentError
    @test (
        occursin("variable_order", sprint(showerror, err)) ||
        occursin("at least one", sprint(showerror, err))
    )
end
