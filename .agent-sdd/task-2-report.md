# Task 2 Report: Q-Degree Helper And Replay Branch

## Status

Done.

## TDD Evidence

### RED

Command:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Observed before implementation:

- Exit code: `1`
- Failure: `UndefVarError: normalize_sl3_local_q_degree not defined in Suslin`
- Location: `test/expert/sl3_local_q_degree_normalization.jl:66`

This matched the expected missing implementation state from Task 1.

### GREEN

Command:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Observed after implementation:

- Exit code: `0`
- Summary: `local SL3 q-degree normalization certificate | 51 passed`

## Implementation Summary

### `src/algorithm/sl3_local.jl`

- Factored the existing recognizer into:
  - `_recognize_sl3_local_special_form_matrix`
  - `_recognize_sl3_local_special_form_parameters`
  - `_recognize_sl3_local_supported_family`
- Kept `realize_sl3_local` and `realize_sl3_local_certificate` behavior intact by routing them through the new special-form recognizer and existing family selection.
- Added:
  - `normalize_sl3_local_q_degree(A, X; check_monic=true)`
  - `normalize_sl3_local_q_degree(p, q, r, s, X; check_monic=true)`
  - `_normalize_sl3_local_q_degree_form(form)`
- Used polynomial `divrem(form.q, form.p)` for quotient/remainder construction.
- Added ordinary-polynomial and variable-degree helpers:
  - `_sl3_local_require_ordinary_polynomial_ring`
  - `_sl3_local_variable_index`
  - `_sl3_local_degree_in_variable`
- Extended certificate verification for `branch == :q_degree_normalization` with:
  - partial replay reconstruction
  - witness-key validation
  - witness algebra checks
  - expected-factor matching against the single correction factor

### `test/expert/sl3_local_q_degree_normalization.jl`

- Kept scope limited to the controller-requested alignment fix:
  - made the helper call near line 133 pass `check_monic = true` explicitly

## Tests Run

1. `julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'`
   - Before implementation: failed with missing method
   - After implementation: passed
2. `julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'`
   - Passed

## Files Changed

- `src/algorithm/sl3_local.jl`
- `test/expert/sl3_local_q_degree_normalization.jl`
- `.agent-sdd/task-2-report.md`

## Self-Review

- The new normalization path is isolated to a dedicated certificate branch and does not alter the factor construction for the existing supported families.
- Verifier changes are narrowly scoped: normal branches still use `verify_factorization`, while q-degree normalization uses reconstruction against the normalized target and correction factor.
- The implementation follows the brief’s requirement to use Oscar/AbstractAlgebra polynomial division via `divrem`.
- No public driver, registration, or unrelated files were modified.
