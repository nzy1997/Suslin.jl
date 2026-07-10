using Suslin
using Test

include("ci/TestManifest.jl")
include("ci/TestRunner.jl")
using .TestManifest
using .TestRunner

const TEST_MANIFEST = load_manifest(joinpath(@__DIR__, "ci", "shards.toml"))
validate_manifest(TEST_MANIFEST, @__DIR__)

for (target, files) in requested_targets(copy(ARGS), TEST_MANIFEST)
    @testset "$target" begin
        for file in files
            TestRunner.timed_test_file(file) do
                include(file)
            end
        end
    end
end
