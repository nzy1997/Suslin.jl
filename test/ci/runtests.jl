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

function load_error_message(toml::AbstractString)
    return mktemp() do path, io
        write(io, toml)
        flush(io)
        try
            load_manifest(path)
        catch error
            return sprint(showerror, error)
        end
        return nothing
    end
end

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
        raw"public\api_surface.jl",
        raw"\\server\share\api_surface.jl",
        "//server/share/api_surface.jl",
    ]
    for path in invalid_paths
        @test_throws ArgumentError TestManifest.validate_relative_path(path, "test path")
    end
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
