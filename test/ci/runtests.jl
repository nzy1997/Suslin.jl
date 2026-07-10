using Test

include("TestManifest.jl")
using .TestManifest

const TEST_ROOT = normpath(joinpath(@__DIR__, ".."))
const MANIFEST_PATH = joinpath(@__DIR__, "shards.toml")
const EXPECTED_SHARDS = [
    "public",
    "internal-core",
    "internal-fixtures",
    "expert-core",
    "expert-laurent-a",
    "expert-laurent-b",
    "expert-sl3",
    "expert-quillen",
    "expert-ecp",
    "expert-integration",
]

@testset "CI shard manifest" begin
    manifest = load_manifest(MANIFEST_PATH)
    @test validate_manifest(manifest, TEST_ROOT) === nothing
    @test shard_ids(manifest) == EXPECTED_SHARDS
    @test length(files_for_group(manifest, "public")) == 6
    @test length(files_for_group(manifest, "internal")) == 31
    @test length(files_for_group(manifest, "expert")) == 85
    @test length(all_test_files(manifest)) == 122
    @test length(unique(all_test_files(manifest))) == 122
    @test owner_shard(manifest, "expert/documentation_smoke.jl") == "expert-integration"
    @test manifest.documentation_smoke == "expert/documentation_smoke.jl"
end
