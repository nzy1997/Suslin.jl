# Issue 155 Lazy Laurent Peel No Initial Determinant Design

## Context

Issue #155 adds a guardrail for the lazy determinant route. The current
Laurent column peel validates determinant-one input by calling
`classify_laurent_determinant(A)` before the first peel step. The Laurent
`GL_n` normalization path also classifies the original input determinant up
front. Later performance work needs a behavioral test proving a lazy entry
point can peel first and defer determinant classification until a smaller
submatrix is available.

Issue #154 is present on `main`: it added
`test/fixtures/laurent_lazy_determinant_cases.jl` and the internal fixture
validator. The smallest useful supported fixture is
`monomial-unit-row-column-cores`, a `3 x 3` Laurent matrix with a monomial-unit
determinant.

## Approaches

Recommended approach: add a focused internal lazy peel surface,
`Suslin._factor_laurent_gl_lazy_determinant_peel(A; determinant_probe,
progress_callback)`, and a focused expert test that injects a sentinel
determinant probe. The lazy surface validates only shape and ring, performs one
real column peel, emits a completed progress event, then calls the determinant
probe on the smaller block. If the deferred determinant is not `:one`, it throws
a staged-boundary error. This proves the ordering without implementing the full
lazy correction certificate.

Alternative 1: add only the test against a not-yet-implemented function. This
would satisfy the "red" side but would not deliver the required green command.

Alternative 2: retrofit `normalize_laurent_gl_matrix(A)` with a probe hook.
This would exercise the eager normalization API instead of the future lazy peel
entry point and would not isolate the behavior the issue asks for.

## Design

Create `test/expert/laurent_lazy_peel_no_initial_det.jl`. The test loads the
#154 lazy determinant fixtures, selects the monomial-unit fixture, and records
progress callback events. Its lazy sentinel fails with a recognizable error if
it is called at the original matrix dimension or before `completed_steps >= 1`.
The test expects the lazy route to call the sentinel only after a completed peel
and on a smaller block.

Add a negative eager control in the same file. The control calls an internal
eager shim via `_factor_laurent_sl_column_peel(A; determinant_probe = ...)`.
That sentinel raises a recognizable error on the original matrix. The test
asserts the eager path probes at the original dimension with no completed peel
step, proving the guard would fail under the old eager determinant behavior.

Implementation changes stay internal to `src/algorithm/laurent_column_peel.jl`.
The existing eager peel accepts an optional `determinant_probe` keyword whose
default is `classify_laurent_determinant`, preserving existing callers. The new
lazy entry point shares the existing peel-step and certificate construction
helpers. It does not export a public API and does not implement determinant
correction for non-`:one` deferred profiles; it throws a staged-boundary
`ArgumentError` after demonstrating the deferred probe order.

Register the expert test in `test/runtests.jl` so full expert runs include the
guard.

## Testing

Use TDD:

- First add `test/expert/laurent_lazy_peel_no_initial_det.jl` and register it.
- Run
  `julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'`
  and confirm it fails because the lazy entry point and eager probe keyword are
  not implemented.
- Implement the internal lazy/eager probe surfaces.
- Re-run the focused command, then run
  `julia --project=. -e 'using Pkg; Pkg.test()'`.

The test must assert that the first determinant probe observed by the lazy path
sees a matrix smaller than the original and that at least one completed peel
step was recorded before the probe.
