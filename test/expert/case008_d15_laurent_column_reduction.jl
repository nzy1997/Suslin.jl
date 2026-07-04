using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl"))

function _case008_d15_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case008_d15_apply_factors(factors, column, R)
    n = length(column)
    return _case008_d15_reduction_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _case008_d15_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _case008_d15_tamper_first_factor(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _case008_d15_tamper_certificate_first_factor(cert)
    tampered = _case008_d15_tamper_first_factor(
        cert.factors,
        cert.ring,
        length(cert.original_column),
    )
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        cert.stages,
        tampered,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d15_tamper_stage_coefficient(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stage = stages[stage_idx]
    coefficients = collect(stage.coefficients)
    coefficients[1] += one(cert.ring)
    stages[stage_idx] = merge(stage, (; coefficient = coefficients[1], coefficients = tuple(coefficients...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d15_tamper_stage_source(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stage = stages[stage_idx]
    source_indices = collect(stage.source_indices)
    source_indices[1] = source_indices[1] == 2 ? 3 : 2
    stages[stage_idx] = merge(stage, (; source_index = source_indices[1], source_indices = tuple(source_indices...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d15_tamper_stage_source_type(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stage = stages[stage_idx]
    source_indices = collect(stage.source_indices)
    source_indices[1] = Float64(source_indices[1])
    stages[stage_idx] = merge(
        stage,
        (; source_index = Float64(stage.source_index), source_indices = tuple(source_indices...)),
    )
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d15_tamper_stage_metadata(cert, replacement)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stages[stage_idx] = merge(stages[stage_idx], replacement)
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d15_diagnostic_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

@testset "case_008 d=15 Laurent column reduction" begin
    fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    target = _case008_d15_target_column(R, length(column))

    factors = Suslin.reduce_unimodular_column(column, R)
    @test _case008_d15_apply_factors(factors, column, R) == target

    tampered_factors = _case008_d15_tamper_first_factor(factors, R, length(column))
    @test _case008_d15_apply_factors(tampered_factors, column, R) != target

    certificate = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(certificate)
    @test any(stage -> stage.kind == :laurent_elementary_row_preconditioning, certificate.stages)
    @test !any(stage -> stage.kind == :case008_special_case, certificate.stages)
    @test _case008_d15_apply_factors(
        certificate.factors,
        certificate.original_column,
        certificate.ring,
    ) == target
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_certificate_first_factor(certificate),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_coefficient(certificate),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_source(certificate),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_source_type(certificate),
    )

    preconditioning_stage = certificate.stages[end]
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_metadata(certificate, (; source_index = 2.0)),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_metadata(
            certificate,
            (; source_indices = collect(preconditioning_stage.source_indices)),
        ),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_metadata(certificate, (; source_indices = ())),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_metadata(
            certificate,
            (; coefficients = (preconditioning_stage.coefficients[1], one(R))),
        ),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_metadata(certificate, (; source_index = 3)),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_metadata(
            certificate,
            (; coefficient = preconditioning_stage.coefficient + one(R)),
        ),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_metadata(
            certificate,
            (; source_index = 1, source_indices = (1,)),
        ),
    )

    @test preconditioning_stage.kind == :laurent_elementary_row_preconditioning
    @test preconditioning_stage.target_index == 1
    @test preconditioning_stage.source_indices == Tuple(2:15)
    @test length(preconditioning_stage.coefficients) == 14
    @test preconditioning_stage.coefficient_strategy == :target_unit_laurent_linear_synthesis
    @test preconditioning_stage.transformed_certificate.stages[end].kind == :unit_entry

    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 15
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test :laurent_elementary_row_preconditioning in diagnostic.attempted_stages
    @test !(:case008_special_case in diagnostic.attempted_stages)
    detail = _case008_d15_diagnostic_stage_detail(
        diagnostic,
        :laurent_elementary_row_preconditioning,
    )
    @test detail !== nothing
    @test detail.outcome == :supported
    @test detail.target_index == 1
    @test detail.source_indices == Tuple(2:15)
    @test detail.coefficient_strategy == :target_unit_laurent_linear_synthesis
    @test detail.coefficient_count == 14
    @test detail.transformed_stage == :unit_entry

    negative = ToricBuilderCase008D15ColumnBoundary.non_unimodular_negative_control(fixture)
    negative_diagnostic = Suslin.diagnose_unimodular_column_reduction(
        negative.failing_column,
        negative.ring,
    )
    @test negative_diagnostic.status == :precondition_failed
    @test negative_diagnostic.failure_code == :not_unimodular
    @test isempty(negative_diagnostic.attempted_stages)
    @test_throws ArgumentError Suslin.reduce_unimodular_column(
        negative.failing_column,
        negative.ring,
    )
end

@testset "case_008 d=15 Laurent row-preconditioning guard coverage" begin
    R, (u, _) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["guard_u", "guard_v"])
    column = [idx == 1 ? zero(R) : idx == 2 ? one(R) : zero(R) for idx in 1:15]
    column16 = [idx == 1 ? zero(R) : idx == 10 ? one(R) : zero(R) for idx in 1:16]

    @test Suslin._laurent_row_preconditioning_synthesis_coefficients(
        column,
        R,
        1,
        (),
    ) === nothing
    @test Suslin._laurent_row_preconditioning_synthesis_coefficients(
        column,
        R,
        1,
        (2,);
        solver = (A, B) -> error("No exact solution exists for A * U = B"),
    ) === nothing
    @test_throws ErrorException Suslin._laurent_row_preconditioning_synthesis_coefficients(
        column,
        R,
        1,
        (2,);
        solver = (A, B) -> error("unexpected row-preconditioning solve failure"),
    )

    fixed_mismatch = (;
        target_index = 1,
        source_indices = (2, 3),
        coefficient_strategy = :fixed_coefficients,
        coefficients = (one(R),),
        max_nonzero_coefficients = 2,
    )
    unknown_strategy = (;
        target_index = 1,
        source_indices = (2,),
        coefficient_strategy = :unknown_strategy,
        coefficients = (one(R),),
        max_nonzero_coefficients = 1,
    )
    @test Suslin._laurent_row_preconditioning_coefficients(
        column,
        R,
        fixed_mismatch,
    ) === nothing
    @test Suslin._laurent_row_preconditioning_coefficients(
        column,
        R,
        unknown_strategy,
    ) === nothing

    @test !Suslin._laurent_row_preconditioning_source_order_ok((), Tuple(2:15))
    @test !Suslin._laurent_row_preconditioning_source_order_ok((3, 2), Tuple(2:15))
    @test Suslin._laurent_row_preconditioning_source_order_ok((2, 4), Tuple(2:15))

    @test Suslin._laurent_row_preconditioning_target_unit_equation_ok(
        column,
        R,
        1,
        (2,),
        (one(R),),
    )
    @test !Suslin._laurent_row_preconditioning_target_unit_equation_ok(
        column,
        R,
        1,
        (2,),
        (u,),
    )

    @test Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        (2,),
        (one(R),),
        :target_unit_laurent_linear_synthesis,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        (2,),
        (one(R), one(R)),
        :target_unit_laurent_linear_synthesis,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        (),
        (),
        :target_unit_laurent_linear_synthesis,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        (2,),
        (zero(R),),
        :target_unit_laurent_linear_synthesis,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        2,
        (3,),
        (one(R),),
        :target_unit_laurent_linear_synthesis,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        (2,),
        (one(R),),
        :fixed_coefficients,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        Tuple(2:16),
        ntuple(_ -> one(R), 15),
        :target_unit_laurent_linear_synthesis,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        (3, 2),
        (one(R), one(R)),
        :target_unit_laurent_linear_synthesis,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column,
        R,
        1,
        (2,),
        (u,),
        :target_unit_laurent_linear_synthesis,
    )

    @test Suslin._laurent_row_preconditioning_stage_spec_ok(
        column16,
        R,
        1,
        (10,),
        (one(R),),
        :fixed_coefficients,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column16,
        R,
        1,
        (9,),
        (one(R),),
        :fixed_coefficients,
    )
    @test !Suslin._laurent_row_preconditioning_stage_spec_ok(
        column16,
        R,
        1,
        (10,),
        (u,),
        :fixed_coefficients,
    )
end
