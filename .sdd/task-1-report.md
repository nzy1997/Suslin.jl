# Task 1 Report: Failing Validator for Issue 38 Fixture

## What I implemented
- Added `test/internal/toricbuilder_issue38_fixture.jl`.
- Defined the issue-38 fixture path constant and the required-field tuple.
- Implemented internal validation helpers for:
  - metadata and required fields
  - determinant profile checks
  - row normalization checks
  - column normalization checks
  - expected current-status / failure-message metadata
- Added `validate_toricbuilder_issue38_fixture(entry)` and `validate_toricbuilder_issue38_catalog(catalog)`.
- Added the required testset that first checks for the fixture file, then includes it and validates the catalog.

## Tests run and results
- `julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'`
  - First attempt failed because the Julia environment was not instantiated and `Oscar` was unavailable.
- `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
  - Succeeded and installed/precompiled the project dependencies locally.
- `julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'`
  - Failed as expected because `test/fixtures/toricbuilder_issue38_cases.jl` does not exist yet.

## TDD evidence
- RED command:
  - `julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'`
- Relevant failure:
  - `Expression: isfile(TORICBUILDER_ISSUE38_FIXTURE_PATH)`
  - `SystemError: opening file ".../test/fixtures/toricbuilder_issue38_cases.jl": No such file or directory`

## Files changed
- `test/internal/toricbuilder_issue38_fixture.jl`
- `.sdd/task-1-report.md`

## Self-review findings
- The validator is scoped to the planned future fixture shape and does not add any fixture data.
- The testset is intentionally red right now because the fixture file is missing.
- The current `expected_current_status` check is permissive about the exact inner shape so task 2 can still choose a concrete representation.

## Concerns
- The repo needed `Pkg.instantiate()` before the exact Julia command could reach the missing-fixture failure.
- The validator assumes the future fixture will expose the metadata and normalization objects described in the task briefs.
