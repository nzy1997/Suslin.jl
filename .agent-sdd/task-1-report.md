# Task 1 Report: Focused Q-Degree Normalization Test

## Scope and Result
- Created `test/expert/sl3_local_q_degree_normalization.jl` exactly for the focused q-degree normalization certificate behavior described in `task-1-brief.md`, with only syntax-level preservation for Julia execution.
- Kept edits limited to test ownership files:
  - `test/expert/sl3_local_q_degree_normalization.jl`
  - `.agent-sdd/task-1-report.md` (report artifact)
- Did not modify source implementation, registration, runners, or docs.

## Focused Execution
- Command run:
  - `julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'`
- Expected RED confirmed.
- Failure observed:
  - `UndefVarError: normalize_sl3_local_q_degree not defined in Suslin`
  - Trigger location: `_assert_q_degree_certificate` call in the testset.

## Self-review of test quality
- Real behavior checks:
  - Verifies certificate type, branch, target, selected variable, witness fields, witness degree bookkeeping (`p_degree`, `remainder_degree`), and relation `q = p * quotient + remainder`.
  - Verifies normalization identity via `normalized_target * correction == target`.
  - Verifies both scalar-entry API and matrix-entry API call paths (`normalize_sl3_local_q_degree(...)` both tuple form and matrix form).
  - Verifies certificate replay with `Suslin.verify_sl3_local_realization`.
- Fixture reuse:
  - Uses `SL3MurthyGuptaFixtureCatalog.catalog()` for the canonical `mg-q-degree-normalization` case and adds an independent hand-checked case.
- Determinant-one negative controls:
  - Includes non-monic leading coefficient check (`@test_throws ArgumentError`) and corruption tests that perturb witness `remainder` and add an unexpected witness field; both assert verifier failure.
- Public-driver isolation:
  - No runner registration, no public test suite wiring changes, and no implementation files changed.
  - Only this expert test file and report were added.
