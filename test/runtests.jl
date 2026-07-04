using Suslin
using Test

const TEST_GROUP_FILES = Dict(
    "public" => [
        "public/api_surface.jl",
        "public/laurent_gl_certificate_options.jl",
        "public/factorization_driver_shell.jl",
        "public/park_woodburn_polynomial_factorization.jl",
        "public/laurent_large_acceptance.jl",
        "public/toricbuilder_factor_toric_block_acceptance.jl",
    ],
    "internal" => [
        "internal/rings.jl",
        "internal/laurent_rings.jl",
        "internal/laurent_fixtures.jl",
        "internal/laurent_normalization.jl",
        "internal/laurent_linear_solve.jl",
        "internal/gl_laurent_normalization.jl",
        "internal/laurent_lazy_determinant_fixtures.jl",
        "internal/toricbuilder_contract.jl",
        "internal/toricbuilder_issue38_fixture.jl",
        "internal/toricbuilder_cache_q_blocks.jl",
        "internal/toricbuilder_case010_column_boundary.jl",
        "internal/toricbuilder_case008_d21_column_boundary.jl",
        "internal/toricbuilder_case008_d16_column_boundary.jl",
        "internal/toricbuilder_case008_d15_column_boundary.jl",
        "internal/toricbuilder_cache_status_report.jl",
        "internal/toricbuilder_cache_case010_certificate.jl",
        "internal/toricbuilder_problem_catalog.jl",
        "internal/sl3_murthy_gupta_fixtures.jl",
        "internal/ecp_column_fixtures.jl",
        "internal/ecp_mainline_fixtures.jl",
        "internal/quillen_patch_fixtures.jl",
        "internal/quillen_mainline_fixtures.jl",
        "internal/park_woodburn_polynomial_fixtures.jl",
        "internal/park_woodburn_sl3_driver_fixtures.jl",
        "internal/park_woodburn_sln_driver_fixtures.jl",
        "internal/park_woodburn_mainline_acceptance_fixtures.jl",
        "internal/steinberg_optimization_fixtures.jl",
    ],
    "expert" => [
        "expert/elementary_matrices.jl",
        "expert/steinberg_factor_count_optimization.jl",
        "expert/elementary_preconditioning.jl",
        "expert/block_embeddings.jl",
        "expert/documentation_smoke.jl",
        "expert/factorization_small_examples.jl",
        "expert/cohn_type.jl",
        "expert/normality.jl",
        "expert/normality_rank_one.jl",
        "expert/polynomial_normality_fixtures.jl",
        "expert/polynomial_normality_support_boundary.jl",
        "expert/laurent_elementary_core.jl",
        "expert/sl3_local.jl",
        "expert/sl3_local_extended.jl",
        "expert/sl3_local_certificate.jl",
        "expert/sl3_local_local_factors.jl",
        "expert/sl3_local_q_degree_normalization.jl",
        "expert/sl3_local_split_lemma.jl",
        "expert/sl3_local_murthy_q_unit.jl",
        "expert/sl3_local_murthy_resultant.jl",
        "expert/sl3_local_murthy_context.jl",
        "expert/sl3_local_murthy_gupta.jl",
        "expert/sln_to_sl3_reduction.jl",
        "expert/park_woodburn_route_certificate.jl",
        "expert/park_woodburn_sl3_driver_context.jl",
        "expert/park_woodburn_sl3_witness_selection.jl",
        "expert/park_woodburn_sl3_local_evidence_provider.jl",
        "expert/park_woodburn_sln_peel_step.jl",
        "expert/park_woodburn_polynomial_column_peel.jl",
        "expert/park_woodburn_sln_recursive_driver.jl",
        "expert/park_woodburn_sln_driver_context.jl",
        "expert/sln_to_sl3_diagnostics.jl",
        "expert/laurent_column_peel_issue38.jl",
        "expert/laurent_lazy_peel_no_initial_det.jl",
        "expert/laurent_lazy_peel_certificate.jl",
        "expert/laurent_lazy_submatrix_normalization.jl",
        "expert/laurent_lazy_correction_hoist.jl",
        "expert/laurent_lazy_row_column_correction.jl",
        "expert/laurent_column_reduction_diagnostics.jl",
        "expert/case008_d21_laurent_column_reduction.jl",
        "expert/case008_d15_laurent_column_reduction.jl",
        "expert/case008_d16_laurent_column_reduction.jl",
        "expert/case010_laurent_column_reduction.jl",
        "expert/issue38_laurent_gl_certificate.jl",
        "expert/quillen_induction.jl",
        "expert/quillen_patching_exact.jl",
        "expert/quillen_denominator_cover.jl",
        "expert/quillen_denominator_cover_solver.jl",
        "expert/quillen_patch_substitution_chain.jl",
        "expert/quillen_supplied_evidence_patch_assembly.jl",
        "expert/quillen_local_factor_sequence.jl",
        "expert/quillen_denominator_extraction.jl",
        "expert/quillen_local_certificate.jl",
        "expert/quillen_contribution_normalization.jl",
        "expert/quillen_global_patch_assembly.jl",
        "expert/quillen_patch_verification_hardening.jl",
        "expert/quillen_induction_constructive.jl",
        "expert/park_woodburn_quillen_route_adapter.jl",
        "expert/quillen_murthy_adapter_consumption.jl",
        "expert/unimodular_columns.jl",
        "expert/unimodular_reduction_exact.jl",
        "expert/ecp_input_context.jl",
        "expert/ecp_column_certificate.jl",
        "expert/ecp_variable_change_replay.jl",
        "expert/ecp_monicity_normalization.jl",
        "expert/ecp_monicity_search.jl",
        "expert/ecp_link_witness_general.jl",
        "expert/ecp_link_witnesses.jl",
        "expert/ecp_link_step.jl",
        "expert/ecp_link_step_general.jl",
        "expert/ecp_induction_normality_general.jl",
        "expert/ecp_induction_normality.jl",
        "expert/elementary_column_property.jl",
    ],
)

function requested_test_groups(args::Vector{String})
    if isempty(args)
        return ["public", "internal"]
    end

    groups = String[]
    for arg in args
        append!(groups, filter(!isempty, split(arg, ',')))
    end

    if "all" in groups
        return ["public", "internal", "expert"]
    end

    invalid = sort(setdiff(unique(groups), collect(keys(TEST_GROUP_FILES))))
    isempty(invalid) || throw(ArgumentError("unknown test groups: $(join(invalid, ", "))"))

    return unique(groups)
end

for group in requested_test_groups(copy(ARGS))
    @testset "$group" begin
        for file in TEST_GROUP_FILES[group]
            include(file)
        end
    end
end
