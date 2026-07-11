# Issue 355 Task 1 Report

## Status

Task 1 is complete. This change adds the focused red contract for the case_008
d14 Laurent-to-polynomial bridge and registers it in the CI test manifest. No
production or report implementation logic was added.

## Changes

- Added `test/expert/case008_d14_laurent_to_polynomial_bridge.jl`.
- Registered the test as `expert-core` in `test/ci/shards.toml` and the
  manifest expectation list in `test/ci/runtests.jl`.
- Updated manifest count assertions from 126 to 127 total tests and from 86 to
  87 expert tests.
- Left `test/runtests.jl` unchanged because it loads test paths through the CI
  manifest rather than maintaining a separate explicit list.

The contract validates the d14 fixture, checks registration, requires the
specified report fields, checks the expected case metadata and successful final
status, and calls the future report validator.

## Verification

`git diff --check` passed.

`julia --project=. test/ci/runtests.jl` passed: 15 testsets, 188 tests.

`julia --project=. -e 'include("test/expert/case008_d14_laurent_to_polynomial_bridge.jl")'`
failed as intended with:

`UndefVarError: case008_d14_laurent_to_polynomial_bridge_report not defined`

The focused run passed the fixture and registration checks before reaching the
missing future report constructor. The failure is therefore the required red
state, not a syntax or fixture-loading failure.

## Scope Notes

The pre-existing modification to
`test/internal/laurent_to_polynomial_certificate.jl` and the pre-existing plan
file were preserved and are not part of this task's implementation.
