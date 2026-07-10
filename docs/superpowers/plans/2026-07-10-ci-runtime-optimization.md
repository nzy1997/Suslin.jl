# CI Runtime Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 111-minute per-PR covered full suite with fail-closed affected-test selection, fresh 99% patch coverage, and a complete coverage run every day at 06:00 Asia/Shanghai.

**Architecture:** A checked-in manifest assigns every currently registered test to one of ten stable shards and maps source/fixture changes to affected shards. Pull requests run at most four selected shards and aggregate their LCOV reports into one `PR Gate`; a separate scheduled workflow runs all ten shards and refreshes Codecov's carried-forward `full-suite` baseline.

**Tech Stack:** Julia 1.12, Julia `Test` and `TOML` stdlibs, GitHub Actions matrices, `julia-actions`, LCOV, `diff-cover` 10.1.0, Codecov flags.

## Global Constraints

- Preserve the existing `test/runtests.jl` interfaces: no arguments, `public`, `internal`, `expert`, and `all` keep their current meanings.
- Use exactly these coverage shard ids: `public`, `internal-core`, `internal-fixtures`, `expert-core`, `expert-laurent-a`, `expert-laurent-b`, `expert-sl3`, `expert-quillen`, `expert-ecp`, and `expert-integration`.
- A source-changing PR must never produce an empty target set; selection errors and unknown source paths select all ten shards.
- Pull-request and complete-coverage matrices use `max-parallel: 4`.
- Changed executable source lines require at least 99% fresh patch coverage.
- The complete workflow runs at `0 22 * * *` UTC, supports manual dispatch, and runs on tag pushes.
- The complete workflow contains every file currently registered by `test/runtests.jl all`, with no omissions or duplicates.
- Codecov `full-suite` uses carryforward; `pr-selected` never uses carryforward.
- Keep the current Codecov project `target: auto` and 1% threshold.
- Do not change GitHub branch protection automatically; document `PR Gate` as the required-check follow-up.
- Do not add a package dependency to `Project.toml`; test tooling uses Julia stdlibs plus the external `diff-cover` executable.
- Preserve unrelated staged and untracked files in the working tree. Every commit in this plan stages only the files named by its task.

---

## File Structure

### New files

- `test/ci/shards.toml` — ordered test catalog, shard assignments, complete-run triggers, documentation paths, source impacts, and fixture impacts.
- `test/ci/TestManifest.jl` — load, validate, and query the manifest without running package tests.
- `test/ci/TestRunner.jl` — translate legacy group arguments and new target arguments into ordered file lists.
- `test/ci/TestSelection.jl` — pure changed-path-to-target selection plus JSON encoding for Actions matrices.
- `test/ci/select_shards.jl` — CLI that obtains a Git diff and emits either line output or GitHub outputs.
- `test/ci/check_patch_coverage.sh` — shared local/CI `diff-cover` wrapper with the fixed 99% threshold.
- `test/ci/coverage_changed.sh` — local changed-test and temporary-LCOV driver.
- `test/ci/runtests.jl` — fast unit and negative-control tests for the manifest, runner, selector, and matrix encoding.
- `test/ci/check_patch_coverage_test.sh` — synthetic Git/LCOV positive and negative controls.
- `.github/workflows/Nightly.yml` — scheduled, manual, and tag-triggered complete coverage matrix.

### Modified files

- `test/runtests.jl` — load the manifest, preserve legacy groups, add `shard:<id>` and `documentation-smoke`, and print per-file timings.
- `.github/workflows/CI.yml` — replace duplicate PR jobs with selection, selected shard coverage, and aggregate `PR Gate`; preserve manual documentation behavior.
- `codecov.yml` — declare carryforward behavior for `full-suite` and `pr-selected`.
- `README.md` — document legacy commands, local changed coverage, PR behavior, complete behavior, schedule, and freshness limitation.

---

### Task 1: Add the ordered shard manifest and validator

**Files:**
- Create: `test/ci/shards.toml`
- Create: `test/ci/TestManifest.jl`
- Create: `test/ci/runtests.jl`

**Interfaces:**
- Produces: `TestManifest.Manifest`, `load_manifest(path)::Manifest`, `validate_manifest(manifest, test_root)::Nothing`, `shard_ids(manifest)::Vector{String}`, `files_for_shard(manifest, id)::Vector{String}`, `files_for_group(manifest, group)::Vector{String}`, `all_test_files(manifest)::Vector{String}`, and `owner_shard(manifest, test_path)::Union{Nothing,String}`.
- Consumes: Julia `TOML` and `Test` stdlibs only.

- [ ] **Step 1: Write failing manifest tests**

Create `test/ci/runtests.jl` with this first testset:

```julia
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
```

- [ ] **Step 2: Run the tests and verify the missing module fails**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
```

Expected: `LoadError` because `test/ci/TestManifest.jl` does not exist.

- [ ] **Step 3: Create the manifest loader and validator**

Create `test/ci/TestManifest.jl` with these exact public types and functions:

```julia
module TestManifest

using TOML

export Manifest
export all_test_files
export files_for_group
export files_for_shard
export load_manifest
export owner_shard
export shard_ids
export validate_manifest

const VALID_GROUPS = Set(["public", "internal", "expert"])

struct TestEntry
    path::String
    group::String
    shard::String
end

struct Manifest
    shard_order::Vector{String}
    tests::Vector{TestEntry}
    documentation_smoke::String
    full_run_paths::Set{String}
    full_run_prefixes::Vector{String}
    documentation_paths::Set{String}
    documentation_prefixes::Vector{String}
    source_impacts::Dict{String,Vector{String}}
    fixture_impacts::Dict{String,Vector{String}}
end

function string_vector(value, label::AbstractString)
    value isa Vector || throw(ArgumentError("$label must be an array"))
    all(item -> item isa AbstractString, value) ||
        throw(ArgumentError("$label must contain only strings"))
    return String[String(item) for item in value]
end

function string_map(value, label::AbstractString)
    value isa AbstractDict || throw(ArgumentError("$label must be a table"))
    result = Dict{String,Vector{String}}()
    for (key, entries) in value
        result[String(key)] = string_vector(entries, "$label.$key")
    end
    return result
end

function load_manifest(path::AbstractString)
    raw = TOML.parsefile(path)
    get(raw, "version", nothing) == 1 ||
        throw(ArgumentError("unsupported shard manifest version"))

    shard_order = string_vector(raw["shard_order"], "shard_order")
    tests = TestEntry[]
    for item in raw["tests"]
        push!(tests, TestEntry(
            String(item["path"]),
            String(item["group"]),
            String(item["shard"]),
        ))
    end

    return Manifest(
        shard_order,
        tests,
        String(raw["documentation_smoke"]),
        Set(string_vector(raw["full_run_paths"], "full_run_paths")),
        string_vector(raw["full_run_prefixes"], "full_run_prefixes"),
        Set(string_vector(raw["documentation_paths"], "documentation_paths")),
        string_vector(raw["documentation_prefixes"], "documentation_prefixes"),
        string_map(raw["source_impacts"], "source_impacts"),
        string_map(raw["fixture_impacts"], "fixture_impacts"),
    )
end

shard_ids(manifest::Manifest) = copy(manifest.shard_order)
all_test_files(manifest::Manifest) = [entry.path for entry in manifest.tests]
files_for_group(manifest::Manifest, group::AbstractString) =
    [entry.path for entry in manifest.tests if entry.group == group]
files_for_shard(manifest::Manifest, shard::AbstractString) =
    [entry.path for entry in manifest.tests if entry.shard == shard]

function owner_shard(manifest::Manifest, test_path::AbstractString)
    normalized = startswith(test_path, "test/") ? test_path[6:end] : String(test_path)
    index = findfirst(entry -> entry.path == normalized, manifest.tests)
    return isnothing(index) ? nothing : manifest.tests[index].shard
end

function validate_manifest(manifest::Manifest, test_root::AbstractString)
    length(unique(manifest.shard_order)) == length(manifest.shard_order) ||
        throw(ArgumentError("shard ids must be unique"))
    paths = all_test_files(manifest)
    length(unique(paths)) == length(paths) ||
        throw(ArgumentError("test paths must be unique"))

    shard_set = Set(manifest.shard_order)
    for entry in manifest.tests
        entry.group in VALID_GROUPS ||
            throw(ArgumentError("invalid test group: $(entry.group)"))
        entry.shard in shard_set ||
            throw(ArgumentError("unknown shard $(entry.shard) for $(entry.path)"))
        isfile(joinpath(test_root, entry.path)) ||
            throw(ArgumentError("missing test file: $(entry.path)"))
    end

    for shard in manifest.shard_order
        isempty(files_for_shard(manifest, shard)) &&
            throw(ArgumentError("empty shard: $shard"))
    end

    manifest.documentation_smoke in paths ||
        throw(ArgumentError("documentation smoke test must belong to the complete suite"))

    for impacts in values(manifest.source_impacts)
        all(shard -> shard in shard_set, impacts) ||
            throw(ArgumentError("source impact references an unknown shard"))
    end
    for impacts in values(manifest.fixture_impacts)
        all(shard -> shard in shard_set, impacts) ||
            throw(ArgumentError("fixture impact references an unknown shard"))
    end
    return nothing
end

end
```

- [ ] **Step 4: Create the complete ordered manifest**

Create `test/ci/shards.toml` with `version = 1`, the ten ids from
`EXPECTED_SHARDS` in `shard_order`, and
`documentation_smoke = "expert/documentation_smoke.jl"`.

Add one `[[tests]]` table for every path currently listed in
`TEST_GROUP_FILES`, in the same public/internal/expert order. Each table has
`path`, `group`, and `shard`. Use these exact assignments:

- `public`: `public/api_surface.jl`,
  `public/laurent_gl_certificate_options.jl`,
  `public/factorization_driver_shell.jl`,
  `public/park_woodburn_polynomial_factorization.jl`,
  `public/laurent_large_acceptance.jl`, and
  `public/toricbuilder_factor_toric_block_acceptance.jl`.
- `internal-core`: `internal/rings.jl`, `internal/laurent_rings.jl`,
  `internal/laurent_normalization.jl`, `internal/laurent_linear_solve.jl`,
  `internal/gl_laurent_normalization.jl`, `internal/toricbuilder_contract.jl`,
  `internal/laurent_descent_measure_helpers.jl`,
  `internal/laurent_link_witness_helpers.jl`,
  `internal/laurent_endpoint_reduction_helpers.jl`,
  `internal/toricbuilder_cache_status_report.jl`,
  `internal/toricbuilder_cache_case010_certificate.jl`,
  `internal/sl3_murthy_gupta_fixtures.jl`,
  `internal/ecp_column_fixtures.jl`, `internal/quillen_patch_fixtures.jl`, and
  `internal/park_woodburn_sl3_driver_fixtures.jl`.
- `internal-fixtures`: the other sixteen current internal entries:
  `internal/laurent_fixtures.jl`,
  `internal/laurent_lazy_determinant_fixtures.jl`,
  `internal/toricbuilder_issue38_fixture.jl`,
  `internal/toricbuilder_cache_q_blocks.jl`,
  `internal/toricbuilder_case010_column_boundary.jl`,
  `internal/toricbuilder_case008_d21_column_boundary.jl`,
  `internal/toricbuilder_case008_d16_column_boundary.jl`,
  `internal/toricbuilder_case008_d15_column_boundary.jl`,
  `internal/toricbuilder_case008_d14_column_boundary.jl`,
  `internal/toricbuilder_problem_catalog.jl`,
  `internal/ecp_mainline_fixtures.jl`,
  `internal/quillen_mainline_fixtures.jl`,
  `internal/park_woodburn_polynomial_fixtures.jl`,
  `internal/park_woodburn_sln_driver_fixtures.jl`,
  `internal/park_woodburn_mainline_acceptance_fixtures.jl`, and
  `internal/steinberg_optimization_fixtures.jl`.
- `expert-core`: `expert/elementary_matrices.jl`,
  `expert/steinberg_factor_count_optimization.jl`,
  `expert/elementary_preconditioning.jl`, `expert/block_embeddings.jl`,
  `expert/factorization_small_examples.jl`, `expert/cohn_type.jl`,
  `expert/normality.jl`, `expert/normality_rank_one.jl`,
  `expert/unimodular_columns.jl`, and
  `expert/unimodular_reduction_exact.jl`.
- `expert-integration`: `expert/documentation_smoke.jl`,
  `expert/polynomial_normality_fixtures.jl`,
  `expert/polynomial_normality_support_boundary.jl`,
  `expert/park_woodburn_route_certificate.jl`, and
  `expert/park_woodburn_polynomial_column_peel.jl`.
- `expert-sl3`: every current expert path from `expert/sl3_local.jl` through
  `expert/sln_to_sl3_reduction.jl`; `expert/park_woodburn_sl3_driver_context.jl`,
  `expert/park_woodburn_sl3_witness_selection.jl`,
  `expert/park_woodburn_sl3_local_evidence_provider.jl`,
  `expert/park_woodburn_sln_peel_step.jl`,
  `expert/park_woodburn_sln_recursive_driver.jl`,
  `expert/park_woodburn_sln_driver_context.jl`, and
  `expert/sln_to_sl3_diagnostics.jl`.
- `expert-quillen`: every current expert path from
  `expert/quillen_induction.jl` through
  `expert/quillen_murthy_adapter_consumption.jl`.
- `expert-ecp`: every current expert path from `expert/ecp_input_context.jl`
  through `expert/elementary_column_property.jl`.
- `expert-laurent-a`: `expert/laurent_elementary_core.jl`,
  `expert/laurent_lazy_peel_no_initial_det.jl`,
  `expert/laurent_lazy_submatrix_normalization.jl`,
  `expert/laurent_lazy_row_column_correction.jl`,
  `expert/laurent_native_ecp_boundary_diagnostics.jl`,
  `expert/case008_d14_laurent_descent_measure_contract.jl`,
  `expert/laurent_descent_step_certificate.jl`,
  `expert/case008_d14_laurent_link_witness_context.jl`,
  `expert/laurent_link_witness_certificate.jl`,
  `expert/case008_d14_laurent_endpoint_reduction_search.jl`,
  `expert/case008_d21_laurent_column_reduction.jl`,
  `expert/case008_d16_laurent_column_reduction.jl`, and
  `expert/issue38_laurent_gl_certificate.jl`.
- `expert-laurent-b`: `expert/laurent_column_peel_issue38.jl`,
  `expert/laurent_lazy_peel_certificate.jl`,
  `expert/laurent_lazy_correction_hoist.jl`,
  `expert/laurent_column_reduction_diagnostics.jl`,
  `expert/case008_d14_laurent_descent_profile.jl`,
  `expert/case008_d14_laurent_elementary_move_search.jl`,
  `expert/case008_d14_laurent_post_descent_profile.jl`,
  `expert/case008_d14_laurent_link_witness_search.jl`,
  `expert/case008_d14_laurent_endpoint_reduction_context.jl`,
  `expert/laurent_endpoint_reduction_certificate.jl`,
  `expert/case008_d15_laurent_column_reduction.jl`, and
  `expert/case010_laurent_column_reduction.jl`.

Add these policy values:

```toml
full_run_paths = [
  "Project.toml",
  "Manifest.toml",
  "codecov.yml",
  "src/Suslin.jl",
  "test/runtests.jl",
]
full_run_prefixes = ["src/core/", "test/ci/", ".github/workflows/"]
documentation_paths = ["README.md"]
documentation_prefixes = ["docs/"]
```

Add exact source impacts. Broad, highly shared algorithms intentionally select
all ten shards:

```toml
[source_impacts]
"src/algorithm/cohn_type.jl" = ["public", "expert-core"]
"src/algorithm/column_reduction.jl" = ["public", "internal-core", "internal-fixtures", "expert-core", "expert-laurent-a", "expert-laurent-b", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"]
"src/algorithm/column_reduction_case010.jl" = ["public", "internal-fixtures", "expert-laurent-b"]
"src/algorithm/factorization.jl" = ["public", "internal-core", "internal-fixtures", "expert-core", "expert-laurent-a", "expert-laurent-b", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"]
"src/algorithm/laurent_column_peel.jl" = ["public", "internal-core", "internal-fixtures", "expert-laurent-a", "expert-laurent-b", "expert-integration"]
"src/algorithm/laurent_gl_certificate.jl" = ["public", "internal-core", "internal-fixtures", "expert-laurent-a", "expert-laurent-b"]
"src/algorithm/normality.jl" = ["public", "expert-core", "expert-ecp", "expert-integration"]
"src/algorithm/polynomial_column_peel.jl" = ["public", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"]
"src/algorithm/quillen_induction.jl" = ["public", "internal-core", "internal-fixtures", "expert-quillen", "expert-ecp", "expert-integration"]
"src/algorithm/redundancy.jl" = ["public", "expert-core"]
"src/algorithm/sl3_local.jl" = ["public", "internal-core", "internal-fixtures", "expert-sl3", "expert-quillen", "expert-integration"]
"src/algorithm/sln_to_sl3_reduction.jl" = ["public", "expert-sl3", "expert-integration"]
```

Map each current `test/fixtures/*.jl` path to all internal/expert shards that
consume its domain. Use these exact conservative families:

```toml
[fixture_impacts]
"test/fixtures/laurent_cases.jl" = ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"]
"test/fixtures/laurent_large_acceptance_cases.jl" = ["public", "expert-integration"]
"test/fixtures/laurent_lazy_determinant_cases.jl" = ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"]
"test/fixtures/polynomial_normality_cases.jl" = ["expert-core", "expert-integration"]
"test/fixtures/quillen_mainline_cases.jl" = ["internal-fixtures", "expert-quillen", "expert-integration"]
"test/fixtures/quillen_patch_cases.jl" = ["internal-core", "expert-quillen"]
"test/fixtures/ecp_column_cases.jl" = ["internal-core", "expert-ecp"]
"test/fixtures/ecp_mainline_cases.jl" = ["internal-fixtures", "expert-ecp", "expert-integration"]
"test/fixtures/sl3_murthy_gupta_cases.jl" = ["internal-core", "expert-sl3", "expert-quillen"]
"test/fixtures/park_woodburn_polynomial_cases.jl" = ["internal-fixtures", "expert-quillen", "expert-ecp", "expert-integration"]
"test/fixtures/park_woodburn_sl3_driver_cases.jl" = ["internal-core", "expert-sl3", "expert-quillen", "expert-integration"]
"test/fixtures/park_woodburn_sln_driver_cases.jl" = ["internal-fixtures", "expert-sl3", "expert-ecp", "expert-integration"]
"test/fixtures/park_woodburn_mainline_acceptance_cases.jl" = ["internal-fixtures", "expert-sl3", "expert-quillen", "expert-ecp", "expert-integration"]
"test/fixtures/steinberg_optimization_cases.jl" = ["internal-fixtures", "expert-core"]
"test/fixtures/toricbuilder_cache_q_blocks.jl" = ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"]
"test/fixtures/toricbuilder_case008_d14_column_boundary.jl" = ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"]
"test/fixtures/toricbuilder_case008_d15_column_boundary.jl" = ["internal-fixtures", "expert-laurent-b"]
"test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl" = ["expert-laurent-b"]
"test/fixtures/toricbuilder_case008_d16_column_boundary.jl" = ["internal-fixtures", "expert-laurent-a"]
"test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl" = ["expert-laurent-a"]
"test/fixtures/toricbuilder_case008_d21_column_boundary.jl" = ["internal-fixtures", "expert-laurent-a"]
"test/fixtures/toricbuilder_case010_column_boundary.jl" = ["internal-fixtures", "expert-laurent-b"]
"test/fixtures/toricbuilder_factor_toric_block_3.jl" = ["public", "expert-integration"]
"test/fixtures/toricbuilder_issue38_cases.jl" = ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"]
"test/fixtures/toricbuilder_laurent_problem_catalog.jl" = ["internal-fixtures", "expert-laurent-a", "expert-laurent-b"]
```

- [ ] **Step 5: Run manifest tests**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
```

Expected: `CI shard manifest` passes with 122 unique tests across ten non-empty
shards.

- [ ] **Step 6: Commit the manifest component**

```bash
git add test/ci/shards.toml test/ci/TestManifest.jl test/ci/runtests.jl
git commit -m "test: define CI shard manifest"
```

---

### Task 2: Make the test runner shard-aware and timed

**Files:**
- Create: `test/ci/TestRunner.jl`
- Modify: `test/ci/runtests.jl`
- Modify: `test/runtests.jl:1-160`

**Interfaces:**
- Consumes: `TestManifest.Manifest`, `files_for_group`, `files_for_shard`.
- Produces: `TestRunner.requested_targets(args, manifest)::Vector{Pair{String,Vector{String}}}` and the CLI targets `shard:<id>` and `documentation-smoke`.

- [ ] **Step 1: Add failing runner selection tests**

Append to `test/ci/runtests.jl`:

```julia
include("TestRunner.jl")
using .TestRunner

@testset "CI test target resolution" begin
    manifest = load_manifest(MANIFEST_PATH)
    default_targets = requested_targets(String[], manifest)
    @test first.(default_targets) == ["public", "internal"]
    @test length(last(default_targets[1])) == 6
    @test length(last(default_targets[2])) == 31

    all_targets = requested_targets(["all"], manifest)
    @test first.(all_targets) == ["public", "internal", "expert"]
    @test sum(length(last(target)) for target in all_targets) == 122

    shard_target = only(requested_targets(["shard:expert-quillen"], manifest))
    @test first(shard_target) == "shard:expert-quillen"
    @test last(shard_target) == files_for_shard(manifest, "expert-quillen")

    smoke_target = only(requested_targets(["documentation-smoke"], manifest))
    @test last(smoke_target) == [manifest.documentation_smoke]
    @test_throws ArgumentError requested_targets(["shard:missing"], manifest)
    @test_throws ArgumentError requested_targets(["not-a-group"], manifest)
end
```

- [ ] **Step 2: Run and verify the missing runner module fails**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
```

Expected: `LoadError` for `test/ci/TestRunner.jl`.

- [ ] **Step 3: Implement pure target resolution**

Create `test/ci/TestRunner.jl`:

```julia
module TestRunner

using ..TestManifest

export requested_targets

function requested_targets(args::Vector{String}, manifest::Manifest)
    names = isempty(args) ? ["public", "internal"] :
        reduce(vcat, [filter(!isempty, split(arg, ',')) for arg in args])
    "all" in names && (names = ["public", "internal", "expert"])

    targets = Pair{String,Vector{String}}[]
    seen = Set{String}()
    for name in names
        name in seen && continue
        push!(seen, name)
        if name in ("public", "internal", "expert")
            push!(targets, name => files_for_group(manifest, name))
        elseif name == "documentation-smoke"
            push!(targets, name => [manifest.documentation_smoke])
        elseif startswith(name, "shard:")
            shard = name[7:end]
            shard in shard_ids(manifest) ||
                throw(ArgumentError("unknown test shard: $shard"))
            push!(targets, name => files_for_shard(manifest, shard))
        else
            throw(ArgumentError("unknown test target: $name"))
        end
    end
    return targets
end

end
```

- [ ] **Step 4: Replace the hard-coded runner catalog with the manifest**

In `test/runtests.jl`, retain `using Suslin` and `using Test`, then load the two
CI modules and replace `TEST_GROUP_FILES`, `requested_test_groups`, and the
final loop with:

```julia
include("ci/TestManifest.jl")
include("ci/TestRunner.jl")
using .TestManifest
using .TestRunner

const TEST_MANIFEST = load_manifest(joinpath(@__DIR__, "ci", "shards.toml"))
validate_manifest(TEST_MANIFEST, @__DIR__)

for (target, files) in requested_targets(copy(ARGS), TEST_MANIFEST)
    @testset "$target" begin
        for file in files
            elapsed = @elapsed include(file)
            println("TEST_FILE_TIME\t", file, "\t", round(elapsed; digits=3))
        end
    end
end
```

- [ ] **Step 5: Verify legacy and new interfaces**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
julia --startup-file=no --project=. test/runtests.jl documentation-smoke
julia --startup-file=no --project=. test/runtests.jl shard:public
```

Expected: CI unit tests pass; documentation smoke passes; public shard reports
six `TEST_FILE_TIME` lines and passes 729 tests or the current equivalent.

- [ ] **Step 6: Commit the shard-aware runner**

```bash
git add test/ci/TestRunner.jl test/ci/runtests.jl test/runtests.jl
git commit -m "test: add timed shard targets"
```

---

### Task 3: Implement fail-closed affected-test selection

**Files:**
- Create: `test/ci/TestSelection.jl`
- Create: `test/ci/select_shards.jl`
- Modify: `test/ci/runtests.jl`

**Interfaces:**
- Consumes: `Manifest`, `owner_shard`, `shard_ids`.
- Produces: `Selection(targets::Vector{String}, documentation_only::Bool, reasons::Vector{String})`, `select_targets(changed_paths, manifest)::Selection`, `matrix_json(targets)::String`, and CLI options `--base`, `--head`, `--format=lines`, and `--github-output`.

- [ ] **Step 1: Write selector positive and negative controls**

Append to `test/ci/runtests.jl`:

```julia
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
```

- [ ] **Step 2: Run and verify selection symbols are missing**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
```

Expected: `LoadError` for `test/ci/TestSelection.jl`.

- [ ] **Step 3: Implement the pure selector**

Create `test/ci/TestSelection.jl`:

```julia
module TestSelection

using ..TestManifest

export Selection
export matrix_json
export select_targets

struct Selection
    targets::Vector{String}
    documentation_only::Bool
    reasons::Vector{String}
end

is_prefixed(path::AbstractString, prefixes) =
    any(prefix -> startswith(path, prefix), prefixes)

function ordered_targets(selected::Set{String}, manifest::Manifest)
    return [shard for shard in shard_ids(manifest) if shard in selected]
end

function select_targets(changed_paths::Vector{String}, manifest::Manifest)
    isempty(changed_paths) && return Selection(
        shard_ids(manifest), false, ["empty diff: full fallback"])

    if all(path -> path in manifest.documentation_paths ||
                   is_prefixed(path, manifest.documentation_prefixes), changed_paths)
        return Selection(["documentation-smoke"], true, ["documentation-only diff"])
    end

    selected = Set{String}()
    reasons = String[]
    for path in changed_paths
        if path in manifest.full_run_paths || is_prefixed(path, manifest.full_run_prefixes)
            return Selection(shard_ids(manifest), false, ["full-run trigger: $path"])
        elseif startswith(path, "src/")
            impacts = get(manifest.source_impacts, path, nothing)
            isnothing(impacts) && return Selection(
                shard_ids(manifest), false, ["unknown source path: $path"])
            union!(selected, impacts)
            push!(reasons, "$path => $(join(impacts, ','))")
        elseif startswith(path, "test/fixtures/")
            impacts = get(manifest.fixture_impacts, path, nothing)
            isnothing(impacts) && return Selection(
                shard_ids(manifest), false, ["unknown fixture path: $path"])
            union!(selected, impacts)
            push!(reasons, "$path => $(join(impacts, ','))")
        elseif startswith(path, "test/")
            owner = owner_shard(manifest, path)
            isnothing(owner) && return Selection(
                shard_ids(manifest), false, ["unknown test path: $path"])
            push!(selected, owner)
            push!(reasons, "$path => $owner")
        elseif path in manifest.documentation_paths ||
               is_prefixed(path, manifest.documentation_prefixes)
            push!(reasons, "$path => documentation companion")
        else
            return Selection(shard_ids(manifest), false, ["unknown path: $path"])
        end
    end

    isempty(selected) && return Selection(
        shard_ids(manifest), false, ["empty source selection: full fallback"])
    return Selection(ordered_targets(selected, manifest), false, reasons)
end

function matrix_json(targets::Vector{String})
    escaped = replace.(targets, "\\" => "\\\\", "\"" => "\\\"")
    return "[" * join(["\"$target\"" for target in escaped], ",") * "]"
end

end
```

- [ ] **Step 4: Implement the selector CLI**

Create `test/ci/select_shards.jl`. It must:

1. load `test/ci/shards.toml` and validate it;
2. parse `--base=<rev>`, `--head=<rev>` (default `HEAD`),
   `--format=lines`, and `--github-output=<path>`;
3. obtain paths with
   `git diff --name-only --diff-filter=ACMRT <base>...<head>`;
4. call `select_targets`;
5. print one target per line for line format; and
6. append `matrix=<json>`, `documentation_only=<true|false>`, and a multiline
   `reason` value to the requested GitHub output file.

Use this exact git helper so revision arguments are never evaluated by a shell:

```julia
function changed_paths(base::AbstractString, head::AbstractString)
    range = "$base...$head"
    output = read(`git diff --name-only --diff-filter=ACMRT $range`, String)
    return filter(!isempty, split(chomp(output), '\n'))
end
```

Wrap `main(ARGS)` in `try/catch`, print the complete exception with
`showerror(stderr, err, catch_backtrace())`, and exit nonzero. The workflow,
not the CLI, owns the hard-coded all-shards fallback.

- [ ] **Step 5: Verify pure tests and CLI formats**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
julia --startup-file=no --project=. test/ci/select_shards.jl \
  --base=HEAD~1 --head=HEAD --format=lines
```

Expected: unit tests pass; the CLI prints either selected shard ids or all ten
ids, never an empty result.

- [ ] **Step 6: Commit affected-test selection**

```bash
git add test/ci/TestSelection.jl test/ci/select_shards.jl test/ci/runtests.jl
git commit -m "test: select affected CI shards"
```

---

### Task 4: Add the shared local patch-coverage gate

**Files:**
- Create: `test/ci/check_patch_coverage.sh`
- Create: `test/ci/coverage_changed.sh`
- Create: `test/ci/check_patch_coverage_test.sh`

**Interfaces:**
- Consumes: `select_shards.jl --format=lines`, `test/runtests.jl shard:<id>`, LCOV reports, and `uvx` or `DIFF_COVER_BIN`.
- Produces: `coverage_changed.sh --base=<rev>` and `check_patch_coverage.sh --base=<rev> <report>...`, both returning nonzero below 99%.

- [ ] **Step 1: Write a synthetic failing coverage-gate test**

Create `test/ci/check_patch_coverage_test.sh` with this complete temporary-repo
negative control:

```bash
#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

cd "$tmpdir"
git init -q
git config user.name "CI Coverage Test"
git config user.email "ci-coverage@example.invalid"
mkdir -p src
printf 'x = 1\n' > src/example.jl
git add src/example.jl
git commit -qm "base"
printf 'x = 2\n' > src/example.jl

source_file="$tmpdir/src/example.jl"
printf 'TN:\nSF:%s\nDA:1,1\nend_of_record\n' "$source_file" > covered.info
printf 'TN:\nSF:%s\nDA:1,0\nend_of_record\n' "$source_file" > uncovered.info

DIFF_COVER_BIN="uvx --from diff-cover==10.1.0 diff-cover" \
  "$repo_root/test/ci/check_patch_coverage.sh" --base=HEAD~1 covered.info

if DIFF_COVER_BIN="uvx --from diff-cover==10.1.0 diff-cover" \
  "$repo_root/test/ci/check_patch_coverage.sh" --base=HEAD~1 uncovered.info; then
  echo "expected uncovered patch to fail" >&2
  exit 1
fi
```

- [ ] **Step 2: Run and verify the missing wrapper fails**

Run:

```bash
bash test/ci/check_patch_coverage_test.sh
```

Expected: failure because `test/ci/check_patch_coverage.sh` does not exist.

- [ ] **Step 3: Implement the reusable patch gate**

Create executable `test/ci/check_patch_coverage.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

base="origin/main"
if [[ "${1:-}" == --base=* ]]; then
  base="${1#--base=}"
  shift
fi

if [[ "$#" -eq 0 ]]; then
  echo "no LCOV reports supplied" >&2
  exit 2
fi

if [[ -n "${DIFF_COVER_BIN:-}" ]]; then
  read -r -a diff_cover_command <<<"$DIFF_COVER_BIN"
elif command -v uvx >/dev/null 2>&1; then
  diff_cover_command=(uvx --from diff-cover==10.1.0 diff-cover)
else
  echo "uvx is required; install uv from https://docs.astral.sh/uv/" >&2
  exit 2
fi

"${diff_cover_command[@]}" "$@" \
  --compare-branch="$base" \
  --include='src/**' \
  --show-uncovered \
  --fail-under=99
```

- [ ] **Step 4: Implement the local changed-coverage driver**

Create executable `test/ci/coverage_changed.sh` with this complete driver:

```bash
#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

base="origin/main"
if [[ "${1:-}" == --base=* ]]; then
  base="${1#--base=}"
  shift
fi
if [[ "$#" -ne 0 ]]; then
  echo "usage: $0 [--base=<revision>]" >&2
  exit 2
fi
command -v uvx >/dev/null 2>&1 || {
  echo "uvx is required; install uv from https://docs.astral.sh/uv/" >&2
  exit 2
}

targets="$(${JULIA:-julia} --startup-file=no --project=. \
  test/ci/select_shards.jl --base="$base" --head=HEAD --format=lines)"
[[ -n "$targets" ]] || {
  echo "selector returned no targets" >&2
  exit 2
}

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
reports=()

while IFS= read -r target; do
  [[ -n "$target" ]] || continue
  if [[ "$target" == documentation-smoke ]]; then
    "${JULIA:-julia}" --startup-file=no --project=. \
      test/runtests.jl documentation-smoke
    continue
  fi
  report="$tmpdir/${target}.info"
  "${JULIA:-julia}" --startup-file=no --project=. \
    --code-coverage="$report" test/runtests.jl "shard:$target"
  reports+=("$report")
done <<<"$targets"

if [[ "${#reports[@]}" -eq 0 ]]; then
  exit 0
fi
test/ci/check_patch_coverage.sh --base="$base" "${reports[@]}"
```

- [ ] **Step 5: Run positive and negative controls**

Run:

```bash
chmod +x test/ci/check_patch_coverage.sh \
  test/ci/coverage_changed.sh \
  test/ci/check_patch_coverage_test.sh
bash test/ci/check_patch_coverage_test.sh
```

Expected: the covered synthetic change passes, the uncovered synthetic change
is observed failing by the outer test, and the test script exits zero.

- [ ] **Step 6: Commit local coverage tooling**

```bash
git add test/ci/check_patch_coverage.sh \
  test/ci/coverage_changed.sh \
  test/ci/check_patch_coverage_test.sh
git commit -m "test: add changed-line coverage gate"
```

---

### Task 5: Replace pull-request CI with selected coverage and `PR Gate`

**Files:**
- Modify: `.github/workflows/CI.yml:1-111`

**Interfaces:**
- Consumes: selector line/GitHub output, `shard:<id>`, `documentation-smoke`, and `check_patch_coverage.sh`.
- Produces: matrix check names `Selected Tests (<target>)` and one fixed `PR Gate` result.

- [ ] **Step 1: Add a workflow contract test before editing YAML**

Append a `workflow contract` testset to `test/ci/runtests.jl` that reads
`.github/workflows/CI.yml` as text and asserts it contains:

```julia
workflow = read(joinpath(TEST_ROOT, "..", ".github", "workflows", "CI.yml"), String)
@test occursin("name: PR Gate", workflow)
@test occursin("max-parallel: 4", workflow)
@test occursin("pr-selected", workflow)
@test occursin("test/ci/select_shards.jl", workflow)
@test !occursin("name: Full Suite Tests", workflow)
@test !occursin("name: Default Fast Tests", workflow)
@test !occursin("name: Instantiate Dependencies", workflow)
```

- [ ] **Step 2: Run and verify the workflow contract fails**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
```

Expected: failures for missing `PR Gate`, `max-parallel`, and selector text.

- [ ] **Step 3: Rewrite the pull-request jobs**

Keep the existing workflow name, PR concurrency cancellation, manual
`run_docs` input, and documentation job. Replace `instantiate`, `default-tests`,
and `full-suite-tests` with these jobs:

1. `select-tests` on `ubuntu-latest`, PR-only, checkout with `fetch-depth: 0`,
   setup Julia 1.12, instantiate, validate the manifest, and select targets.
2. `selected-tests`, PR-only, `needs: select-tests`, matrix from
   `fromJSON(needs.select-tests.outputs.matrix)`, `fail-fast: false`, and
   `max-parallel: 4`.
3. `pr-gate`, named exactly `PR Gate`, PR-only, `needs: [select-tests,
   selected-tests]`, and `if: always()`.

The selector step must use a shell-owned fallback that does not depend on Julia
remaining executable:

```bash
all='["public","internal-core","internal-fixtures","expert-core","expert-laurent-a","expert-laurent-b","expert-sl3","expert-quillen","expert-ecp","expert-integration"]'
if julia --startup-file=no --project=. test/ci/select_shards.jl \
    --base='${{ github.event.pull_request.base.sha }}' \
    --head='${{ github.event.pull_request.head.sha }}' \
    --github-output="$GITHUB_OUTPUT"; then
  echo 'selector_fallback=false' >> "$GITHUB_OUTPUT"
else
  echo "matrix=$all" >> "$GITHUB_OUTPUT"
  echo 'documentation_only=false' >> "$GITHUB_OUTPUT"
  echo 'selector_fallback=true' >> "$GITHUB_OUTPUT"
fi
```

For coverage targets, selected test jobs run:

```bash
julia --startup-file=no --project=. --code-coverage=@src \
  test/runtests.jl 'shard:${{ matrix.target }}'
```

Then run `julia-actions/julia-processcoverage@v1`, rename `lcov.info` to
`lcov-${{ matrix.target }}.info`, upload it as artifact
`coverage-${{ matrix.target }}` with one-day retention, and upload the same file
through `codecov/codecov-action@v6` using flag `pr-selected`, a unique name, and
`fail_ci_if_error: false`.

For `documentation-smoke`, run
`julia --startup-file=no --project=. test/runtests.jl documentation-smoke` and
skip coverage processing and uploads.

The aggregate `PR Gate` must first fail when either required job result is not
`success`. For source-changing selections it downloads all `coverage-*`
artifacts into `coverage/`, installs `diff-cover==10.1.0` with
`python -m pip install --disable-pip-version-check diff-cover==10.1.0`, and
runs:

```bash
DIFF_COVER_BIN=diff-cover test/ci/check_patch_coverage.sh \
  --base='${{ github.event.pull_request.base.sha }}' coverage/*.info
```

For documentation-only selections, skip artifact download and patch coverage
after confirming the selected test job succeeded.

- [ ] **Step 4: Validate the workflow contract and YAML syntax**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
ruby -e 'require "yaml"; YAML.load_file(".github/workflows/CI.yml", aliases: true); puts "CI YAML OK"'
```

Expected: CI unit tests pass and Ruby prints `CI YAML OK`.

- [ ] **Step 5: Commit PR CI replacement**

```bash
git add .github/workflows/CI.yml test/ci/runtests.jl
git commit -m "ci: run affected coverage shards on pull requests"
```

---

### Task 6: Add complete nightly coverage and Codecov carryforward

**Files:**
- Create: `.github/workflows/Nightly.yml`
- Modify: `codecov.yml:1-17`
- Modify: `test/ci/runtests.jl`

**Interfaces:**
- Consumes: all ten `shard:<id>` targets.
- Produces: `Complete Coverage (<target>)` matrix jobs and Codecov `full-suite` reports.

- [ ] **Step 1: Add failing nightly and Codecov contract tests**

Append to `test/ci/runtests.jl`:

```julia
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
```

- [ ] **Step 2: Run and verify the missing nightly contract fails**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
```

Expected: nightly file and Codecov flag assertions fail.

- [ ] **Step 3: Create the complete coverage workflow**

Create `.github/workflows/Nightly.yml` with:

- name `Complete Coverage`;
- triggers `schedule` with `cron: '0 22 * * *'`, `workflow_dispatch`, and tag
  pushes matching `'*'`;
- workflow concurrency group `complete-coverage-${{ github.ref }}` with
  `cancel-in-progress: false`;
- one matrix job on `ubuntu-latest`, `fail-fast: false`, `max-parallel: 4`, and
  the exact ten target ids from the manifest;
- checkout v6, Julia setup v3 at 1.12 x64, Julia cache v3, and dependency
  instantiation;
- a 30-minute job timeout;
- `julia --startup-file=no --project=. --code-coverage=@src test/runtests.jl
  "shard:${{ matrix.target }}"`;
- `julia-actions/julia-processcoverage@v1` restricted to `src`; and
- `codecov/codecov-action@v6` with `files: lcov.info`, flag `full-suite`, unique
  name `suslin-full-${{ matrix.target }}`, the existing token, and
  `fail_ci_if_error: true`.

- [ ] **Step 4: Configure explicit Codecov flags**

Append this top-level section to `codecov.yml` without changing the existing
project/patch targets or comment configuration:

```yaml
flags:
  full-suite:
    carryforward: true
  pr-selected:
    carryforward: false
```

- [ ] **Step 5: Validate both workflow files and all contracts**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
ruby -e 'require "yaml"; %w[CI.yml Nightly.yml].each { |f| YAML.load_file(".github/workflows/#{f}", aliases: true) }; puts "workflow YAML OK"'
```

Expected: all CI tests pass and Ruby prints `workflow YAML OK`.

- [ ] **Step 6: Commit complete coverage workflow**

```bash
git add .github/workflows/Nightly.yml codecov.yml test/ci/runtests.jl
git commit -m "ci: refresh complete coverage each morning"
```

---

### Task 7: Document the new contract and run the migration gate

**Files:**
- Modify: `README.md:108-126`
- Test: all files created or modified in Tasks 1–6

**Interfaces:**
- Consumes: all new test and workflow commands.
- Produces: contributor-facing commands and final evidence that the shard union equals the legacy full suite.

- [ ] **Step 1: Update README with exact commands and freshness semantics**

Keep the current four legacy commands and add:

````markdown
For fast coverage feedback on the current branch, run only the affected test
shards and check changed executable lines:

```bash
test/ci/coverage_changed.sh --base=origin/main
```

Pull requests use the same fail-closed selector, run at most four test shards
concurrently, and require at least 99% fresh coverage on changed source lines.
Unknown source changes and shared core changes fall back to every shard.

The complete test and coverage suite runs daily at 06:00 Asia/Shanghai
(`0 22 * * *` UTC), on release tags, and on manual dispatch. Codecov carries
the most recent complete `full-suite` report between scheduled runs, so the
project-wide percentage may be up to 24 hours old; pull-request patch coverage
always comes from the current commit.
````

Add an operational note: after the workflow lands, configure the stable
`PR Gate` check as required branch protection. Do not claim that the repository
configuration performs this automatically.

- [ ] **Step 2: Run fast structural and negative-control verification**

Run:

```bash
julia --startup-file=no --project=. test/ci/runtests.jl
bash test/ci/check_patch_coverage_test.sh
julia --startup-file=no --project=. test/runtests.jl documentation-smoke
ruby -e 'require "yaml"; %w[CI.yml Nightly.yml].each { |f| YAML.load_file(".github/workflows/#{f}", aliases: true) }; YAML.load_file("codecov.yml", aliases: true); puts "all YAML OK"'
git diff --check
```

Expected: all unit/negative controls and documentation smoke pass; Ruby prints
`all YAML OK`; `git diff --check` prints nothing.

- [ ] **Step 3: Run every shard without coverage to prove suite completeness**

Run this exact loop:

```bash
for shard in public internal-core internal-fixtures expert-core \
  expert-laurent-a expert-laurent-b expert-sl3 expert-quillen \
  expert-ecp expert-integration; do
  julia --startup-file=no --project=. test/runtests.jl "shard:$shard" || exit 1
done
```

Expected: all ten shard commands pass. The combined test counts equal the
current `julia --project=. test/runtests.jl all` count; manifest tests already
prove the file-level union is exactly 122 files.

- [ ] **Step 4: Run one local changed-coverage smoke path**

On a branch with a source change, run:

```bash
test/ci/coverage_changed.sh --base=origin/main
```

Expected: selector reasons identify the affected shards, no more than the
mapped targets execute unless a full fallback is justified, and diff-cover
reports at least 99% or identifies exact uncovered changed lines.

- [ ] **Step 5: Commit documentation and final local evidence**

```bash
git add README.md
git commit -m "docs: explain incremental and complete coverage"
```

- [ ] **Step 6: Verify the first hosted baseline before enabling protection**

After the branch is pushed and the workflow is available, manually dispatch
`Complete Coverage` on the branch or its merged commit. Confirm:

1. all ten `Complete Coverage (<target>)` jobs pass;
2. no more than four matrix jobs run simultaneously;
3. Codecov lists all ten uploads under `full-suite` and produces the expected
   project percentage;
4. a source-changing PR produces `pr-selected` uploads and a successful
   `PR Gate`; and
5. a documentation-only PR runs only `documentation-smoke` and still produces
   a successful `PR Gate`.

Only after those five observations should the repository owner configure
`PR Gate` as a required check.

---

## Plan Completion Checks

Before declaring the implementation complete, verify all of the following:

- `git status --short` contains no implementation files outside the intended
  commits and preserves the user's pre-existing plan files.
- `git log --oneline` shows one focused commit for each task component.
- `test/ci/runtests.jl` passes from a fresh Julia process.
- The synthetic uncovered LCOV negative control demonstrably fails inside its
  harness.
- All ten shard commands pass and cover exactly the legacy 122 registered test
  files.
- Both workflow YAML files and `codecov.yml` parse.
- PR jobs use at most four-way parallelism and the complete workflow uses at
  most four-way parallelism.
- The cron expression is exactly `0 22 * * *`.
- `PR Gate` enforces 99% fresh changed-line coverage independently of Codecov
  availability.
- Codecov project coverage retains `target: auto` and the 1% tolerance.
