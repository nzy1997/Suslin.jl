using Test
using TOML

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
const EXPECTED_TESTS = [
    ("public/api_surface.jl", "public", "public"),
    ("public/laurent_gl_certificate_options.jl", "public", "public"),
    ("public/factorization_driver_shell.jl", "public", "public"),
    ("public/park_woodburn_polynomial_factorization.jl", "public", "public"),
    ("public/laurent_large_acceptance.jl", "public", "public"),
    ("public/toricbuilder_factor_toric_block_acceptance.jl", "public", "public"),
    ("internal/rings.jl", "internal", "internal-core"),
    ("internal/laurent_rings.jl", "internal", "internal-core"),
    ("internal/laurent_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/laurent_normalization.jl", "internal", "internal-core"),
    ("internal/laurent_linear_solve.jl", "internal", "internal-core"),
    ("internal/gl_laurent_normalization.jl", "internal", "internal-core"),
    ("internal/laurent_lazy_determinant_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/laurent_to_polynomial_route_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/laurent_noether_certificate.jl", "internal", "internal-core"),
    ("internal/laurent_to_polynomial_certificate.jl", "internal", "internal-core"),
    ("internal/toricbuilder_contract.jl", "internal", "internal-core"),
    ("internal/toricbuilder_issue38_fixture.jl", "internal", "internal-fixtures"),
    ("internal/toricbuilder_cache_q_blocks.jl", "internal", "internal-fixtures"),
    ("internal/toricbuilder_case010_column_boundary.jl", "internal", "internal-fixtures"),
    ("internal/toricbuilder_case008_d21_column_boundary.jl", "internal", "internal-fixtures"),
    ("internal/toricbuilder_case008_d16_column_boundary.jl", "internal", "internal-fixtures"),
    ("internal/toricbuilder_case008_d15_column_boundary.jl", "internal", "internal-fixtures"),
    ("internal/toricbuilder_case008_d14_column_boundary.jl", "internal", "internal-fixtures"),
    ("internal/laurent_descent_measure_helpers.jl", "internal", "internal-core"),
    ("internal/laurent_link_witness_helpers.jl", "internal", "internal-core"),
    ("internal/laurent_endpoint_reduction_helpers.jl", "internal", "internal-core"),
    ("internal/toricbuilder_cache_status_report.jl", "internal", "internal-core"),
    ("internal/toricbuilder_cache_case010_certificate.jl", "internal", "internal-core"),
    ("internal/toricbuilder_problem_catalog.jl", "internal", "internal-fixtures"),
    ("internal/sl3_murthy_gupta_fixtures.jl", "internal", "internal-core"),
    ("internal/ecp_column_fixtures.jl", "internal", "internal-core"),
    ("internal/ecp_mainline_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/quillen_patch_fixtures.jl", "internal", "internal-core"),
    ("internal/quillen_mainline_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/park_woodburn_polynomial_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/park_woodburn_sl3_driver_fixtures.jl", "internal", "internal-core"),
    ("internal/park_woodburn_sln_driver_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/park_woodburn_mainline_acceptance_fixtures.jl", "internal", "internal-fixtures"),
    ("internal/steinberg_optimization_fixtures.jl", "internal", "internal-fixtures"),
    ("expert/elementary_matrices.jl", "expert", "expert-core"),
    ("expert/steinberg_factor_count_optimization.jl", "expert", "expert-core"),
    ("expert/elementary_preconditioning.jl", "expert", "expert-core"),
    ("expert/block_embeddings.jl", "expert", "expert-core"),
    ("expert/documentation_smoke.jl", "expert", "expert-integration"),
    ("expert/factorization_small_examples.jl", "expert", "expert-core"),
    ("expert/cohn_type.jl", "expert", "expert-core"),
    ("expert/normality.jl", "expert", "expert-core"),
    ("expert/normality_rank_one.jl", "expert", "expert-core"),
    ("expert/polynomial_normality_fixtures.jl", "expert", "expert-integration"),
    ("expert/polynomial_normality_support_boundary.jl", "expert", "expert-integration"),
    ("expert/laurent_elementary_core.jl", "expert", "expert-laurent-a"),
    ("expert/sl3_local.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_extended.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_certificate.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_local_factors.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_q_degree_normalization.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_split_lemma.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_murthy_q_unit.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_murthy_resultant.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_murthy_context.jl", "expert", "expert-sl3"),
    ("expert/sl3_local_murthy_gupta.jl", "expert", "expert-sl3"),
    ("expert/sln_to_sl3_reduction.jl", "expert", "expert-sl3"),
    ("expert/park_woodburn_route_certificate.jl", "expert", "expert-integration"),
    ("expert/park_woodburn_sl3_driver_context.jl", "expert", "expert-sl3"),
    ("expert/park_woodburn_sl3_witness_selection.jl", "expert", "expert-sl3"),
    ("expert/park_woodburn_sl3_local_evidence_provider.jl", "expert", "expert-sl3"),
    ("expert/park_woodburn_sln_peel_step.jl", "expert", "expert-sl3"),
    ("expert/park_woodburn_polynomial_column_peel.jl", "expert", "expert-integration"),
    ("expert/park_woodburn_sln_recursive_driver.jl", "expert", "expert-sl3"),
    ("expert/park_woodburn_sln_driver_context.jl", "expert", "expert-sl3"),
    ("expert/sln_to_sl3_diagnostics.jl", "expert", "expert-sl3"),
    ("expert/laurent_column_peel_issue38.jl", "expert", "expert-laurent-b"),
    ("expert/laurent_lazy_peel_no_initial_det.jl", "expert", "expert-laurent-a"),
    ("expert/laurent_lazy_peel_certificate.jl", "expert", "expert-laurent-b"),
    ("expert/laurent_lazy_submatrix_normalization.jl", "expert", "expert-laurent-a"),
    ("expert/laurent_lazy_correction_hoist.jl", "expert", "expert-laurent-b"),
    ("expert/laurent_lazy_row_column_correction.jl", "expert", "expert-laurent-a"),
    ("expert/laurent_column_reduction_diagnostics.jl", "expert", "expert-laurent-b"),
    ("expert/laurent_native_ecp_boundary_diagnostics.jl", "expert", "expert-laurent-a"),
    ("expert/case008_d14_laurent_descent_profile.jl", "expert", "expert-laurent-b"),
    ("expert/case008_d14_laurent_descent_measure_contract.jl", "expert", "expert-laurent-a"),
    ("expert/case008_d14_laurent_elementary_move_search.jl", "expert", "expert-laurent-b"),
    ("expert/laurent_descent_step_certificate.jl", "expert", "expert-laurent-a"),
    ("expert/case008_d14_laurent_post_descent_profile.jl", "expert", "expert-laurent-b"),
    ("expert/case008_d14_laurent_link_witness_context.jl", "expert", "expert-laurent-a"),
    ("expert/case008_d14_laurent_link_witness_search.jl", "expert", "expert-laurent-b"),
    ("expert/laurent_link_witness_certificate.jl", "expert", "expert-laurent-a"),
    ("expert/case008_d14_laurent_endpoint_reduction_context.jl", "expert", "expert-laurent-b"),
    ("expert/case008_d14_laurent_endpoint_reduction_search.jl", "expert", "expert-laurent-a"),
    ("expert/laurent_endpoint_reduction_certificate.jl", "expert", "expert-laurent-b"),
    ("expert/case008_d21_laurent_column_reduction.jl", "expert", "expert-laurent-a"),
    ("expert/case008_d15_laurent_column_reduction.jl", "expert", "expert-laurent-b"),
    ("expert/case008_d16_laurent_column_reduction.jl", "expert", "expert-laurent-a"),
    ("expert/case010_laurent_column_reduction.jl", "expert", "expert-laurent-b"),
    ("expert/issue38_laurent_gl_certificate.jl", "expert", "expert-laurent-a"),
    ("expert/quillen_induction.jl", "expert", "expert-quillen"),
    ("expert/quillen_patching_exact.jl", "expert", "expert-quillen"),
    ("expert/quillen_denominator_cover.jl", "expert", "expert-quillen"),
    ("expert/quillen_denominator_cover_solver.jl", "expert", "expert-quillen"),
    ("expert/quillen_patch_substitution_chain.jl", "expert", "expert-quillen"),
    ("expert/quillen_supplied_evidence_patch_assembly.jl", "expert", "expert-quillen"),
    ("expert/quillen_local_factor_sequence.jl", "expert", "expert-quillen"),
    ("expert/quillen_denominator_extraction.jl", "expert", "expert-quillen"),
    ("expert/quillen_local_certificate.jl", "expert", "expert-quillen"),
    ("expert/quillen_contribution_normalization.jl", "expert", "expert-quillen"),
    ("expert/quillen_global_patch_assembly.jl", "expert", "expert-quillen"),
    ("expert/quillen_patch_verification_hardening.jl", "expert", "expert-quillen"),
    ("expert/quillen_induction_constructive.jl", "expert", "expert-quillen"),
    ("expert/park_woodburn_quillen_route_adapter.jl", "expert", "expert-quillen"),
    ("expert/quillen_murthy_adapter_consumption.jl", "expert", "expert-quillen"),
    ("expert/unimodular_columns.jl", "expert", "expert-core"),
    ("expert/unimodular_reduction_exact.jl", "expert", "expert-core"),
    ("expert/ecp_input_context.jl", "expert", "expert-ecp"),
    ("expert/ecp_column_certificate.jl", "expert", "expert-ecp"),
    ("expert/ecp_variable_change_replay.jl", "expert", "expert-ecp"),
    ("expert/ecp_monicity_normalization.jl", "expert", "expert-ecp"),
    ("expert/ecp_monicity_search.jl", "expert", "expert-ecp"),
    ("expert/ecp_link_witness_general.jl", "expert", "expert-ecp"),
    ("expert/ecp_link_witnesses.jl", "expert", "expert-ecp"),
    ("expert/ecp_link_step.jl", "expert", "expert-ecp"),
    ("expert/ecp_link_step_general.jl", "expert", "expert-ecp"),
    ("expert/ecp_induction_normality_general.jl", "expert", "expert-ecp"),
    ("expert/ecp_induction_normality.jl", "expert", "expert-ecp"),
    ("expert/elementary_column_property.jl", "expert", "expert-ecp"),
    ("expert/laurent_to_polynomial_ecp_bridge.jl", "expert", "expert-ecp"),
]
const EXPECTED_FULL_RUN_PATHS = Set([
    "Manifest.toml",
    "Project.toml",
    "codecov.yml",
    "src/Suslin.jl",
    "test/runtests.jl",
])
const EXPECTED_FULL_RUN_PREFIXES = ["src/core/", "test/ci/", ".github/workflows/"]
const EXPECTED_DOCUMENTATION_PATHS = Set(["README.md"])
const EXPECTED_DOCUMENTATION_PREFIXES = ["docs/"]
const EXPECTED_SOURCE_IMPACTS = Dict(
    "src/algorithm/cohn_type.jl" => ["public", "expert-core"],
    "src/algorithm/column_reduction.jl" => ["public", "internal-core", "internal-fixtures", "expert-core", "expert-laurent-a", "expert-laurent-b", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"],
    "src/algorithm/column_reduction_case010.jl" => ["public", "internal-fixtures", "expert-laurent-b"],
    "src/algorithm/factorization.jl" => ["public", "internal-core", "internal-fixtures", "expert-core", "expert-laurent-a", "expert-laurent-b", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"],
    "src/algorithm/laurent_column_peel.jl" => ["public", "internal-core", "internal-fixtures", "expert-laurent-a", "expert-laurent-b", "expert-integration"],
    "src/algorithm/laurent_gl_certificate.jl" => ["public", "internal-core", "internal-fixtures", "expert-laurent-a", "expert-laurent-b"],
    "src/algorithm/normality.jl" => ["public", "expert-core", "expert-ecp", "expert-integration"],
    "src/algorithm/polynomial_column_peel.jl" => ["public", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"],
    "src/algorithm/quillen_induction.jl" => ["public", "internal-core", "internal-fixtures", "expert-quillen", "expert-ecp", "expert-integration"],
    "src/algorithm/redundancy.jl" => ["public", "expert-core"],
    "src/algorithm/sl3_local.jl" => ["public", "internal-core", "internal-fixtures", "expert-sl3", "expert-quillen", "expert-integration"],
    "src/algorithm/sln_to_sl3_reduction.jl" => ["public", "expert-sl3", "expert-integration"],
)
const EXPECTED_FIXTURE_IMPACTS = Dict(
    "test/fixtures/ecp_column_cases.jl" => ["internal-core", "expert-ecp"],
    "test/fixtures/ecp_mainline_cases.jl" => ["internal-fixtures", "expert-ecp", "expert-integration"],
    "test/fixtures/laurent_cases.jl" => ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"],
    "test/fixtures/laurent_large_acceptance_cases.jl" => ["public", "expert-integration"],
    "test/fixtures/laurent_lazy_determinant_cases.jl" => ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"],
    "test/fixtures/park_woodburn_mainline_acceptance_cases.jl" => ["internal-fixtures", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"],
    "test/fixtures/park_woodburn_polynomial_cases.jl" => ["internal-fixtures", "expert-quillen", "expert-ecp", "expert-integration"],
    "test/fixtures/park_woodburn_sl3_driver_cases.jl" => ["internal-core", "expert-sl3", "expert-quillen", "expert-integration"],
    "test/fixtures/park_woodburn_sln_driver_cases.jl" => ["internal-fixtures", "expert-sl3", "expert-ecp", "expert-integration"],
    "test/fixtures/polynomial_normality_cases.jl" => ["expert-core", "expert-integration"],
    "test/fixtures/quillen_mainline_cases.jl" => ["internal-fixtures", "expert-quillen", "expert-integration"],
    "test/fixtures/quillen_patch_cases.jl" => ["internal-core", "expert-quillen"],
    "test/fixtures/sl3_murthy_gupta_cases.jl" => ["internal-core", "expert-sl3", "expert-quillen"],
    "test/fixtures/steinberg_optimization_cases.jl" => ["internal-fixtures", "expert-core"],
    "test/fixtures/toricbuilder_cache_q_blocks.jl" => ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"],
    "test/fixtures/toricbuilder_case008_d14_column_boundary.jl" => ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"],
    "test/fixtures/toricbuilder_case008_d15_column_boundary.jl" => ["internal-fixtures", "expert-laurent-b"],
    "test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl" => ["expert-laurent-b"],
    "test/fixtures/toricbuilder_case008_d16_column_boundary.jl" => ["internal-fixtures", "expert-laurent-a"],
    "test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl" => ["expert-laurent-a"],
    "test/fixtures/toricbuilder_case008_d21_column_boundary.jl" => ["internal-fixtures", "expert-laurent-a"],
    "test/fixtures/toricbuilder_case010_column_boundary.jl" => ["internal-fixtures", "expert-laurent-b"],
    "test/fixtures/toricbuilder_factor_toric_block_3.jl" => ["public", "expert-integration"],
    "test/fixtures/toricbuilder_issue38_cases.jl" => ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"],
    "test/fixtures/toricbuilder_laurent_problem_catalog.jl" => ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"],
)

function manifest_with(
    manifest::Manifest;
    shard_order = copy(manifest.shard_order),
    tests = copy(manifest.tests),
    documentation_smoke = manifest.documentation_smoke,
    full_run_paths = copy(manifest.full_run_paths),
    full_run_prefixes = copy(manifest.full_run_prefixes),
    documentation_paths = copy(manifest.documentation_paths),
    documentation_prefixes = copy(manifest.documentation_prefixes),
    source_impacts = deepcopy(manifest.source_impacts),
    fixture_impacts = deepcopy(manifest.fixture_impacts),
)
    return Manifest(
        shard_order,
        tests,
        documentation_smoke,
        full_run_paths,
        full_run_prefixes,
        documentation_paths,
        documentation_prefixes,
        source_impacts,
        fixture_impacts,
    )
end

@testset "workflow contract" begin
    workflow = read(joinpath(TEST_ROOT, "..", ".github", "workflows", "CI.yml"), String)
    @test occursin("name: PR Gate", workflow)
    @test occursin("max-parallel: 4", workflow)
    @test occursin("pr-selected", workflow)
    @test occursin("test/ci/select_shards.jl", workflow)
    @test !occursin("name: Full Suite Tests", workflow)
    @test !occursin("name: Default Fast Tests", workflow)
    @test !occursin("name: Instantiate Dependencies", workflow)
end

nightly_path = joinpath(TEST_ROOT, "..", ".github", "workflows", "Nightly.yml")
@testset "complete coverage workflow contract" begin
    @test isfile(nightly_path)
    nightly = isfile(nightly_path) ? read(nightly_path, String) : ""
    @test occursin("cron: '0 22 * * *'", nightly)
    @test occursin("max-parallel: 4", nightly)
    @test occursin("full-suite", nightly)
    @test occursin("workflow_dispatch:", nightly)
    @test occursin("tags: ['*']", nightly)

    codecov = read(joinpath(TEST_ROOT, "..", "codecov.yml"), String)
    @test occursin("full-suite:", codecov)
    @test occursin("carryforward: true", codecov)
    @test occursin("pr-selected:", codecov)
    @test occursin("carryforward: false", codecov)
end

function manifest_with_test_entry(
    manifest::Manifest,
    index::Integer;
    path = manifest.tests[index].path,
    group = manifest.tests[index].group,
    shard = manifest.tests[index].shard,
)
    tests = copy(manifest.tests)
    tests[index] = TestManifest.TestEntry(path, group, shard)
    return manifest_with(manifest; tests)
end

function impacts_with(impacts::Dict{String,Vector{String}}, key, shards)
    result = deepcopy(impacts)
    result[key] = shards
    return result
end

function impacts_without(impacts::Dict{String,Vector{String}}, key)
    result = deepcopy(impacts)
    delete!(result, key)
    return result
end

function test_table_toml(;
    path = "public/api_surface.jl",
    group = "public",
    shard = "public",
)
    fields = (("path", path), ("group", group), ("shard", shard))
    lines = ["[[tests]]"]
    for (name, value) in fields
        isnothing(value) && continue
        rendered = value isa AbstractString ? repr(value) : string(value)
        push!(lines, "$name = $rendered")
    end
    return join(lines, "\n")
end

function minimal_manifest_toml(tests_definition::AbstractString)
    return """
version = 1
shard_order = ["public"]
documentation_smoke = "public/api_surface.jl"
full_run_paths = ["Project.toml"]
full_run_prefixes = ["src/core/"]
documentation_paths = ["README.md"]
documentation_prefixes = ["docs/"]

$tests_definition

[source_impacts]
"src/Suslin.jl" = ["public"]

[fixture_impacts]
"test/fixtures/laurent_cases.jl" = ["public"]
"""
end

function captured_error_message(operation)
    try
        operation()
    catch error
        return sprint(showerror, error)
    end
    return nothing
end

function load_error_message(toml::AbstractString)
    return mktemp() do path, io
        write(io, toml)
        flush(io)
        return captured_error_message(() -> load_manifest(path))
    end
end

@testset "CI shard manifest" begin
    manifest = load_manifest(MANIFEST_PATH)
    @test validate_manifest(manifest, TEST_ROOT) === nothing
    @test shard_ids(manifest) == EXPECTED_SHARDS
    @test length(files_for_group(manifest, "public")) == 6
    @test length(files_for_group(manifest, "internal")) == 34
    @test length(files_for_group(manifest, "expert")) == 86
    @test length(all_test_files(manifest)) == 126
    @test length(unique(all_test_files(manifest))) == 126
    @test owner_shard(manifest, "expert/documentation_smoke.jl") == "expert-integration"
    @test manifest.documentation_smoke == "expert/documentation_smoke.jl"
end

@testset "CI shard manifest exact ownership and policies" begin
    manifest = load_manifest(MANIFEST_PATH)
    expected_paths = [path for (path, _, _) in EXPECTED_TESTS]
    expected_groups = Dict(
        group => [path for (path, entry_group, _) in EXPECTED_TESTS if entry_group == group]
        for group in ("public", "internal", "expert")
    )
    expected_shards = Dict(
        shard => [path for (path, _, entry_shard) in EXPECTED_TESTS if entry_shard == shard]
        for shard in EXPECTED_SHARDS
    )

    @test [(entry.path, entry.group, entry.shard) for entry in manifest.tests] ==
          EXPECTED_TESTS
    @test all_test_files(manifest) == expected_paths
    @test all(files_for_group(manifest, group) == paths for (group, paths) in expected_groups)
    @test all(files_for_shard(manifest, shard) == paths for (shard, paths) in expected_shards)
    @test [owner_shard(manifest, path) for path in expected_paths] ==
          [shard for (_, _, shard) in EXPECTED_TESTS]
    @test manifest.full_run_paths == EXPECTED_FULL_RUN_PATHS
    @test manifest.full_run_prefixes == EXPECTED_FULL_RUN_PREFIXES
    @test manifest.documentation_paths == EXPECTED_DOCUMENTATION_PATHS
    @test manifest.documentation_prefixes == EXPECTED_DOCUMENTATION_PREFIXES
    @test manifest.source_impacts == EXPECTED_SOURCE_IMPACTS
    @test manifest.fixture_impacts == EXPECTED_FIXTURE_IMPACTS
end

@testset "CI shard manifest loader rejects malformed tests" begin
    @test load_error_message(minimal_manifest_toml("")) ==
          "ArgumentError: tests must be an array"
    @test load_error_message(minimal_manifest_toml("tests = \"not an array\"")) ==
          "ArgumentError: tests must be an array"
    @test load_error_message(minimal_manifest_toml("tests = [\"not a table\"]")) ==
          "ArgumentError: tests[1] must be a table"

    valid_table = test_table_toml()
    missing_fields = Dict(
        "path" => test_table_toml(; path = nothing),
        "group" => test_table_toml(; group = nothing),
        "shard" => test_table_toml(; shard = nothing),
    )
    wrong_types = Dict(
        "path" => test_table_toml(; path = 1),
        "group" => test_table_toml(; group = 1),
        "shard" => test_table_toml(; shard = 1),
    )
    for field in ("path", "group", "shard")
        @test load_error_message(
            minimal_manifest_toml("$valid_table\n\n$(missing_fields[field])"),
        ) == "ArgumentError: tests[2].$field is required"
        @test load_error_message(
            minimal_manifest_toml("$valid_table\n\n$(wrong_types[field])"),
        ) == "ArgumentError: tests[2].$field must be a string"
    end

    valid_manifest = minimal_manifest_toml(valid_table)
    @test load_error_message(
        replace(
            valid_manifest,
            "full_run_paths = [\"Project.toml\"]" =>
                "full_run_paths = [\"Project.toml\", \"Project.toml\"]",
        ),
    ) == "ArgumentError: full_run_paths values must be unique"
    @test load_error_message(
        replace(
            valid_manifest,
            "documentation_paths = [\"README.md\"]" =>
                "documentation_paths = [\"README.md\", \"README.md\"]",
        ),
    ) == "ArgumentError: documentation_paths values must be unique"
end

@testset "canonical POSIX repository paths" begin
    @test TestManifest.validate_relative_path(
        "public/api_surface.jl",
        "test path",
    ) === nothing

    invalid_paths = [
        "",
        "public//api_surface.jl",
        "public/./api_surface.jl",
        "public/../internal/rings.jl",
        "/public/api_surface.jl",
        "public/api_surface.jl/",
        "C:/public/api_surface.jl",
        "public/D:/outside.jl",
        "public/D:outside.jl",
        raw"public\api_surface.jl",
        raw"\\server\share\api_surface.jl",
        "//server/share/api_surface.jl",
    ]
    for path in invalid_paths
        @test_throws ArgumentError TestManifest.validate_relative_path(path, "test path")
    end

    manifest = load_manifest(MANIFEST_PATH)
    @test captured_error_message(
        () -> validate_manifest(
            manifest_with(manifest; full_run_paths = Set(["src/D:/outside.jl"])),
            TEST_ROOT,
        ),
    ) == "ArgumentError: full_run_paths entry must not contain ':'"
end

@testset "CI shard manifest validation rejects malformed data" begin
    manifest = load_manifest(MANIFEST_PATH)

    @test_throws ArgumentError validate_manifest(
        manifest_with_test_entry(manifest, 1; path = "public/../internal/rings.jl"),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with_test_entry(manifest, 1; path = "ci/runtests.jl"),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with_test_entry(
            manifest,
            1;
            group = "internal",
            shard = "internal-core",
        ),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with_test_entry(manifest, 1; shard = "internal-core"),
        TEST_ROOT,
    )

    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; full_run_paths = Set{String}()),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; full_run_paths = Set(["missing.file"])),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; full_run_prefixes = String[]),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; full_run_prefixes = ["src/core/", "src/core/"]),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; full_run_prefixes = ["missing/"]),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; documentation_paths = Set{String}()),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; documentation_prefixes = String[]),
        TEST_ROOT,
    )

    @test_throws ArgumentError validate_manifest(
        manifest_with(manifest; source_impacts = Dict{String,Vector{String}}()),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            source_impacts = Dict("test/runtests.jl" => ["public"]),
        ),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            source_impacts = Dict("src/algorithm/not_a_file.jl" => ["public"]),
        ),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            source_impacts = impacts_with(
                manifest.source_impacts,
                "src/algorithm/cohn_type.jl",
                String[],
            ),
        ),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            source_impacts = impacts_with(
                manifest.source_impacts,
                "src/algorithm/cohn_type.jl",
                ["public", "public"],
            ),
        ),
        TEST_ROOT,
    )

    missing_fixture = first(keys(manifest.fixture_impacts))
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            fixture_impacts = impacts_without(manifest.fixture_impacts, missing_fixture),
        ),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            fixture_impacts = impacts_with(
                manifest.fixture_impacts,
                "test/fixtures/not_a_file.jl",
                ["expert-core"],
            ),
        ),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            fixture_impacts = impacts_with(
                manifest.fixture_impacts,
                "test/fixtures/laurent_cases.jl",
                String[],
            ),
        ),
        TEST_ROOT,
    )
    @test_throws ArgumentError validate_manifest(
        manifest_with(
            manifest;
            fixture_impacts = impacts_with(
                manifest.fixture_impacts,
                "test/fixtures/laurent_cases.jl",
                ["internal-fixtures", "internal-fixtures"],
            ),
        ),
        TEST_ROOT,
    )
end

include("TestRunner.jl")
using .TestRunner

@testset "CI test target resolution" begin
    manifest = load_manifest(MANIFEST_PATH)
    default_targets = requested_targets(String[], manifest)
    @test first.(default_targets) == ["public", "internal"]
    @test length(last(default_targets[1])) == 6
    @test length(last(default_targets[2])) == 34

    all_targets = requested_targets(["all"], manifest)
    @test first.(all_targets) == ["public", "internal", "expert"]
    @test sum(length(last(target)) for target in all_targets) == 126
    @test first.(requested_targets(["all,public"], manifest)) ==
          ["public", "internal", "expert"]
    @test_throws ArgumentError requested_targets(["all,not-a-group"], manifest)
    @test_throws ArgumentError requested_targets(["all,shard:missing"], manifest)
    @test_throws ArgumentError requested_targets(["all", "not-a-group"], manifest)
    @test_throws ArgumentError requested_targets(["all", "shard:missing"], manifest)

    shard_target = only(requested_targets(["shard:expert-quillen"], manifest))
    @test first(shard_target) == "shard:expert-quillen"
    @test last(shard_target) == files_for_shard(manifest, "expert-quillen")

    smoke_target = only(requested_targets(["documentation-smoke"], manifest))
    @test last(smoke_target) == [manifest.documentation_smoke]
    @test_throws ArgumentError requested_targets(["shard:missing"], manifest)
    @test_throws ArgumentError requested_targets(["not-a-group"], manifest)
end

@testset "CI test target overlap preflight" begin
    manifest = load_manifest(MANIFEST_PATH)

    duplicate_target = requested_targets(
        ["shard:expert-quillen", "shard:expert-quillen"],
        manifest,
    )
    @test first.(duplicate_target) == ["shard:expert-quillen"]
    @test last(only(duplicate_target)) == files_for_shard(manifest, "expert-quillen")
    @test first.(requested_targets(
        ["documentation-smoke,documentation-smoke"],
        manifest,
    )) == ["documentation-smoke"]

    group_shard_error = captured_error_message(
        () -> requested_targets(["public,shard:public"], manifest),
    )
    @test occursin("public/api_surface.jl", something(group_shard_error, ""))
    @test occursin(
        "both public and shard:public",
        something(group_shard_error, ""),
    )

    smoke_shard_error = captured_error_message(
        () -> requested_targets(
            ["documentation-smoke,shard:expert-integration"],
            manifest,
        ),
    )
    @test occursin(manifest.documentation_smoke, something(smoke_shard_error, ""))
    @test occursin(
        "both documentation-smoke and shard:expert-integration",
        something(smoke_shard_error, ""),
    )

    legacy_targets = requested_targets(["public,internal"], manifest)
    @test first.(legacy_targets) == ["public", "internal"]
    shard_targets = requested_targets(
        ["shard:expert-quillen,shard:expert-ecp"],
        manifest,
    )
    @test first.(shard_targets) == ["shard:expert-quillen", "shard:expert-ecp"]
    @test length(unique(vcat(last.(shard_targets)...))) ==
          sum(length(last(target)) for target in shard_targets)

    @test_throws ArgumentError requested_targets([",", ""], manifest)
    empty_public_manifest = manifest_with(
        manifest;
        tests = filter(entry -> entry.group != "public", manifest.tests),
    )
    @test_throws ArgumentError requested_targets(["public"], empty_public_manifest)
end

timing_lines(output::AbstractString) =
    filter(!isempty, split(chomp(output), '\n'))

@testset "CI test file timing" begin
    normal_io = IOBuffer()
    normal_result = try
        TestRunner.timed_test_file(() -> :completed, "normal.jl"; io = normal_io)
    catch error
        error
    end
    @test normal_result === :completed
    normal_lines = timing_lines(String(take!(normal_io)))
    @test length(normal_lines) == 1
    if length(normal_lines) == 1
        normal_fields = split(only(normal_lines), '\t')
        @test length(normal_fields) == 3
        if length(normal_fields) == 3
            @test normal_fields[1:2] == ["TEST_FILE_TIME", "normal.jl"]
            @test !isnothing(tryparse(Float64, normal_fields[3]))
        end
    end

    test_io = IOBuffer()
    test_result = try
        TestRunner.timed_test_file(() -> (@test true), "test-pass.jl"; io = test_io)
    catch error
        error
    end
    @test test_result isa Test.Pass
    @test length(timing_lines(String(take!(test_io)))) == 1

    expected_error = ErrorException("timed failure")
    error_io = IOBuffer()
    caught_error = try
        TestRunner.timed_test_file(
            () -> throw(expected_error),
            "thrown-error.jl";
            io = error_io,
        )
        nothing
    catch error
        error
    end
    @test caught_error === expected_error
    error_lines = timing_lines(String(take!(error_io)))
    @test length(error_lines) == 1
    if length(error_lines) == 1
        @test startswith(only(error_lines), "TEST_FILE_TIME\tthrown-error.jl\t")
    end

    mktemp() do path, file_io
        write(file_io, "error(\"top-level failure\")\n")
        flush(file_io)
        close(file_io)

        load_error_io = IOBuffer()
        caught_load_error = try
            TestRunner.timed_test_file(
                () -> Base.include(Main, path),
                "temporary-load-error.jl";
                io = load_error_io,
            )
            nothing
        catch error
            error
        end
        @test caught_load_error isa LoadError
        load_error_lines = timing_lines(String(take!(load_error_io)))
        @test length(load_error_lines) == 1
        if length(load_error_lines) == 1
            @test startswith(
                only(load_error_lines),
                "TEST_FILE_TIME\ttemporary-load-error.jl\t",
            )
        end
    end
end

include("TestSelection.jl")
using .TestSelection

@testset "affected test selection" begin
    manifest = load_manifest(MANIFEST_PATH)

    laurent = select_targets(["src/algorithm/laurent_column_peel.jl"], manifest)
    @test laurent.documentation_only == false
    @test laurent.targets == [
        "public", "internal-core", "internal-fixtures",
        "expert-laurent-a", "expert-laurent-b", "expert-integration",
    ]

    sl3 = select_targets(["src/algorithm/sln_to_sl3_reduction.jl"], manifest)
    @test sl3.targets == ["public", "expert-sl3", "expert-integration"]

    owned_test = select_targets(["test/expert/ecp_link_step.jl"], manifest)
    @test owned_test.targets == ["expert-ecp"]

    fixture = select_targets(["test/fixtures/quillen_patch_cases.jl"], manifest)
    @test fixture.targets == ["internal-core", "expert-quillen"]

    docs = select_targets(["README.md", "docs/src/index.md"], manifest)
    @test docs.documentation_only
    @test docs.targets == ["documentation-smoke"]

    core = select_targets(["src/core/rings.jl"], manifest)
    @test core.targets == shard_ids(manifest)

    unknown_source = select_targets(["src/algorithm/new_algorithm.jl"], manifest)
    @test unknown_source.targets == shard_ids(manifest)

    unknown_fixture = select_targets(["test/fixtures/new_fixture.jl"], manifest)
    @test unknown_fixture.targets == shard_ids(manifest)

    empty_diff = select_targets(String[], manifest)
    @test empty_diff.targets == shard_ids(manifest)

    @test matrix_json(["public", "expert-ecp"]) ==
        "[\"public\",\"expert-ecp\"]"
end

@testset "affected test selection ordering and fail-closed controls" begin
    manifest = load_manifest(MANIFEST_PATH)

    combined = select_targets([
        "test/expert/ecp_link_step.jl",
        "src/algorithm/sln_to_sl3_reduction.jl",
        "docs/src/index.md",
    ], manifest)
    @test combined.targets == [
        "public", "expert-sl3", "expert-ecp", "expert-integration",
    ]
    @test !combined.documentation_only
    @test combined.reasons == [
        "test/expert/ecp_link_step.jl => expert-ecp",
        "src/algorithm/sln_to_sl3_reduction.jl => public,expert-sl3,expert-integration",
        "docs/src/index.md => documentation companion",
    ]

    full_run = select_targets(["README.md", "Project.toml"], manifest)
    @test full_run.targets == shard_ids(manifest)
    @test full_run.reasons == ["full-run trigger: Project.toml"]

    unknown_test = select_targets(["test/expert/new_test.jl"], manifest)
    @test unknown_test.targets == shard_ids(manifest)
    @test unknown_test.reasons == ["unknown test path: test/expert/new_test.jl"]

    unknown_path = select_targets(["scripts/new_helper.jl"], manifest)
    @test unknown_path.targets == shard_ids(manifest)
    @test unknown_path.reasons == ["unknown path: scripts/new_helper.jl"]

    @test matrix_json(["a\\b", "c\"d"]) == "[\"a\\\\b\",\"c\\\"d\"]"

    control_targets = [
        "line\nfeed",
        "tab\tvalue",
        "carriage\rreturn",
        "back\bspace",
        "form\ffeed",
        "nul\0value",
        "unit$(Char(0x1f))separator",
        "snow 雪",
    ]
    control_json = matrix_json(control_targets)
    @test control_json ==
        "[\"line\\nfeed\",\"tab\\tvalue\",\"carriage\\rreturn\"," *
        "\"back\\bspace\",\"form\\ffeed\",\"nul\\u0000value\"," *
        "\"unit\\u001fseparator\",\"snow 雪\"]"
    parsed_targets = try
        TOML.parse("targets = $control_json")["targets"]
    catch error
        error
    end
    @test parsed_targets == control_targets
end

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SELECTOR_CLI = joinpath(@__DIR__, "select_shards.jl")

function run_selector_cli(arguments::Vector{String}; directory = REPOSITORY_ROOT)
    stdout = IOBuffer()
    stderr = IOBuffer()
    command = Cmd(
        `$(Base.julia_cmd()) --startup-file=no --project=$REPOSITORY_ROOT $SELECTOR_CLI $arguments`;
        dir = directory,
    )
    process = run(pipeline(ignorestatus(command); stdout, stderr))
    return success(process), String(take!(stdout)), String(take!(stderr))
end

function run_selector_expression(
    expression::AbstractString,
    arguments::Vector{String} = String[];
    directory = REPOSITORY_ROOT,
)
    stdout = IOBuffer()
    stderr = IOBuffer()
    command = Cmd(
        `$(Base.julia_cmd()) --startup-file=no --project=$REPOSITORY_ROOT -e $expression -- $arguments`;
        dir = directory,
    )
    process = run(pipeline(ignorestatus(command); stdout, stderr))
    return success(process), String(take!(stdout)), String(take!(stderr))
end

function changed_paths_for_repository(repository::AbstractString)
    expression = """
    include($(repr(SELECTOR_CLI)))
    foreach(println, changed_paths(
        ARGS[1], ARGS[2]; repository_root = ARGS[3],
    ))
    """
    succeeded, stdout, stderr = run_selector_expression(
        expression,
        ["HEAD~1", "HEAD", String(repository)],
    )
    paths = String.(filter(!isempty, split(chomp(stdout), '\n')))
    return succeeded, paths, stderr
end

function git_command(repository::AbstractString, arguments::Vector{String})
    return run(`git -C $repository $arguments`)
end

function initialize_git_repository(repository::AbstractString)
    git_command(repository, ["init", "--quiet"])
    git_command(repository, ["config", "user.email", "ci-test@example.com"])
    git_command(repository, ["config", "user.name", "CI Test"])
    return nothing
end

function write_repository_file(
    repository::AbstractString,
    path::AbstractString,
    contents::AbstractString,
)
    full_path = joinpath(repository, split(path, '/')...)
    mkpath(dirname(full_path))
    write(full_path, contents)
    return nothing
end

function commit_repository(repository::AbstractString, message::AbstractString)
    git_command(repository, ["add", "--all"])
    git_command(repository, ["commit", "--quiet", "-m", String(message)])
    return nothing
end

@testset "affected test selection CLI" begin
    succeeded, stdout, stderr = run_selector_cli(
        ["--base=HEAD", "--head=HEAD", "--format=lines"];
        directory = joinpath(TEST_ROOT, "expert"),
    )
    @test succeeded
    @test stdout == join(EXPECTED_SHARDS, '\n') * "\n"
    @test isempty(stderr)

    mktemp() do github_output, io
        write(io, "existing=value\n")
        close(io)
        succeeded, stdout, stderr = run_selector_cli([
            "--base=HEAD",
            "--head=HEAD",
            "--github-output=$github_output",
        ])
        @test succeeded
        @test isempty(stdout)
        @test isempty(stderr)
        @test read(github_output, String) ==
            "existing=value\n" *
            "matrix=$(matrix_json(EXPECTED_SHARDS))\n" *
            "documentation_only=false\n" *
            "reason<<EOF\n" *
            "empty diff: full fallback\n" *
            "EOF\n"
    end

    succeeded, stdout, stderr = run_selector_cli(["--format=lines"])
    @test !succeeded
    @test isempty(stdout)
    @test occursin("ArgumentError: --base is required", stderr)
    @test occursin("Stacktrace:", stderr)

    succeeded, stdout, stderr = run_selector_cli([
        "--base=HEAD", "--format=json",
    ])
    @test !succeeded
    @test isempty(stdout)
    @test occursin("ArgumentError: unsupported format: json", stderr)

    succeeded, stdout, stderr = run_selector_cli([
        "--base=HEAD", "--unknown-option",
    ])
    @test !succeeded
    @test isempty(stdout)
    @test occursin("ArgumentError: unknown option: --unknown-option", stderr)
end

@testset "affected test selection GitHub delimiter safety" begin
    mktemp() do github_output, io
        close(io)
        expression = """
        include($(repr(SELECTOR_CLI)))
        selection = Selection(
            ["public"], false, ["first line\nEOF\ninjected=value"],
        )
        write_github_output(ARGS[1], selection)
        """
        succeeded, stdout, stderr = run_selector_expression(
            expression,
            [github_output],
        )
        @test succeeded
        @test isempty(stdout)
        @test isempty(stderr)
        @test read(github_output, String) ==
            "matrix=[\"public\"]\n" *
            "documentation_only=false\n" *
            "reason<<EOF_1\n" *
            "first line\n" *
            "EOF\n" *
            "injected=value\n" *
            "EOF_1\n"
    end
end

@testset "affected test selection Git diff safety" begin
    manifest = load_manifest(MANIFEST_PATH)

    mktempdir() do outside_repository
        succeeded, stdout, stderr = run_selector_cli(
            ["--base=HEAD", "--head=HEAD", "--format=lines"];
            directory = outside_repository,
        )
        @test succeeded
        @test stdout == join(EXPECTED_SHARDS, '\n') * "\n"
        @test isempty(stderr)
    end

    mktempdir() do repository
        initialize_git_repository(repository)
        write_repository_file(repository, "README.md", "initial docs\n")
        write_repository_file(repository, "src/core/rings.jl", "core source\n")
        commit_repository(repository, "initial files")

        write_repository_file(repository, "README.md", "changed docs\n")
        rm(joinpath(repository, "src", "core", "rings.jl"))
        commit_repository(repository, "delete core source")

        succeeded, paths, stderr = changed_paths_for_repository(repository)
        @test succeeded
        @test isempty(stderr)
        @test Set(paths) == Set(["README.md", "src/core/rings.jl"])
        @test select_targets(paths, manifest).targets == shard_ids(manifest)

        write_repository_file(
            repository,
            "src/algorithm/sln_to_sl3_reduction.jl",
            "source rename fixture\n",
        )
        write_repository_file(
            repository,
            "test/expert/ecp_link_step.jl",
            "test rename fixture\n",
        )
        commit_repository(repository, "add rename fixtures")

        mkpath(joinpath(repository, "docs"))
        git_command(repository, [
            "mv", "src/algorithm/sln_to_sl3_reduction.jl", "docs/source_guide.jl",
        ])
        git_command(repository, [
            "mv", "test/expert/ecp_link_step.jl", "docs/test_guide.jl",
        ])
        commit_repository(repository, "rename source and test to docs")

        succeeded, paths, stderr = changed_paths_for_repository(repository)
        @test succeeded
        @test isempty(stderr)
        @test Set(paths) == Set([
            "src/algorithm/sln_to_sl3_reduction.jl",
            "docs/source_guide.jl",
            "test/expert/ecp_link_step.jl",
            "docs/test_guide.jl",
        ])
        @test select_targets(paths, manifest).targets == [
            "public", "expert-sl3", "expert-ecp", "expert-integration",
        ]

        mkpath(joinpath(repository, "src", "algorithm"))
        mkpath(joinpath(repository, "test", "expert"))
        git_command(repository, [
            "mv", "docs/source_guide.jl", "src/algorithm/sln_to_sl3_reduction.jl",
        ])
        git_command(repository, [
            "mv", "docs/test_guide.jl", "test/expert/ecp_link_step.jl",
        ])
        commit_repository(repository, "rename docs to source and test")

        succeeded, paths, stderr = changed_paths_for_repository(repository)
        @test succeeded
        @test isempty(stderr)
        @test Set(paths) == Set([
            "docs/source_guide.jl",
            "src/algorithm/sln_to_sl3_reduction.jl",
            "docs/test_guide.jl",
            "test/expert/ecp_link_step.jl",
        ])

        succeeded, expected_stdout, stderr = run_selector_cli(
            ["--base=HEAD~1", "--head=HEAD", "--format=lines"],
        )
        @test succeeded
        @test !isempty(expected_stdout)
        @test isempty(stderr)

        succeeded, stdout, stderr = run_selector_cli(
            ["--base=HEAD~1", "--head=HEAD", "--format=lines"];
            directory = repository,
        )
        @test succeeded
        @test stdout == expected_stdout
        @test isempty(stderr)
    end
end
