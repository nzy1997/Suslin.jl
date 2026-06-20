using Suslin
using Test

const TEST_GROUP_FILES = Dict(
    "public" => [
        "public/api_surface.jl",
        "public/factorization_driver_shell.jl",
        "public/laurent_large_acceptance.jl",
    ],
    "internal" => [
        "internal/rings.jl",
        "internal/laurent_rings.jl",
        "internal/laurent_fixtures.jl",
        "internal/laurent_normalization.jl",
        "internal/laurent_linear_solve.jl",
        "internal/gl_laurent_normalization.jl",
        "internal/toricbuilder_contract.jl",
        "internal/toricbuilder_issue38_fixture.jl",
    ],
    "expert" => [
        "expert/elementary_matrices.jl",
        "expert/block_embeddings.jl",
        "expert/documentation_smoke.jl",
        "expert/factorization_small_examples.jl",
        "expert/cohn_type.jl",
        "expert/normality.jl",
        "expert/laurent_elementary_core.jl",
        "expert/sl3_local.jl",
        "expert/sl3_local_extended.jl",
        "expert/sln_to_sl3_reduction.jl",
        "expert/quillen_induction.jl",
        "expert/quillen_patching_exact.jl",
        "expert/unimodular_columns.jl",
        "expert/unimodular_reduction_exact.jl",
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
