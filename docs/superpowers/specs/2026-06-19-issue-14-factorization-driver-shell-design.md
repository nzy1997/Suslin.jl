# Issue 14 Factorization Driver Shell Design

## Context

`elementary_factorization(A)` currently mixes public input validation with the
only implemented algorithm slice: a univariate polynomial `3 x 3` local
`SL_3` path. That shape-specific gate blocks valid larger `SL_n` or Laurent
`GL_n` inputs before they can cross the determinant normalization boundary
added for issue #20.

Issues #19 and #6 provide exact Laurent fixture expectations for later
ToricBuilder consumers. Issue #20 exposes `normalize_laurent_gl_matrix(A)`,
which can turn supported Laurent `GL_n` determinants into determinant-one
cores with correction metadata or throw a staged determinant error. Issue #21
defines the test commands: the package test entry point runs default
public/internal tests, and `julia --project=. test/runtests.jl all` is the full
suite.

No visual companion is needed because the change is a public driver control
flow refactor, not a visual design task.

## Design Choice

Refactor `elementary_factorization(A)` into a small driver shell with explicit
layers:

1. Generic matrix validation: require square input and size `n >= 3`.
2. Ring validation: distinguish polynomial rings from Laurent rings, and reject
   unsupported exact-ring categories before algorithm dispatch.
3. Determinant boundary: for Laurent input, call
   `normalize_laurent_gl_matrix(A)` before algorithm support checks; for
   polynomial input, require determinant `1`.
4. Supported-case detection: identify the existing univariate polynomial local
   `SL_3` slice and route it unchanged to `realize_sl3_local`.
5. Staged dispatch failure: valid but unsupported determinant-one inputs fail
   with an algorithm-layer message describing the missing reduction layer.

Recommended approach: keep the current public return type for
`elementary_factorization(A)`. Supported inputs still return an exact factor
sequence. Unsupported normalized Laurent or larger `SL_n` inputs throw staged
`ArgumentError`s instead of returning metadata through the factorization API.
This preserves public API compatibility while still ensuring every Laurent
`GL_n` input crosses the issue #20 boundary first.

Alternative considered: return normalization metadata from
`elementary_factorization(A)` for Laurent `GL_n` inputs. That would satisfy one
reading of the issue text, but it would make one public function return both
factor sequences and non-factor metadata. Later factorization callers would
need new type checks before a stable output contract exists.

Alternative considered: leave Laurent handling in place and only change the
`3 x 3` message. That would not separate generic validation from algorithm
dispatch and would still let nonsquare Laurent inputs fail inside the
normalization layer instead of the generic driver layer.

## Driver Components

The implementation stays in `src/algorithm/factorization.jl` and uses private
helpers:

- `_validate_factorization_matrix(A)`: validates square shape and `n >= 3`.
- `_factorization_ring_profile(R)`: classifies the base ring as
  `:polynomial` or `:laurent`.
- `_normalize_factorization_input(A, profile)`: calls
  `normalize_laurent_gl_matrix(A)` for Laurent rings and returns the matrix
  that should be checked by the `SL_n` algorithm layer.
- `_require_sl_determinant(A)`: rejects polynomial inputs whose determinant is
  not `1`.
- `_is_supported_local_sl3_slice(A, R)`: identifies the current implemented
  univariate local `SL_3` slice.
- `_throw_staged_factorization_failure(A, profile, normalization)`: reports
  the missing algorithm layer for valid determinant-one inputs that are not in
  the current slice.

These helpers are internal; no new public exports are required.

## Error Semantics

Generic matrix errors happen before determinant normalization:

- Nonsquare input throws `ArgumentError("A must be square")`.
- Size `< 3` throws an `ArgumentError` that says the public driver requires
  size at least `3`.

Polynomial determinant errors happen before unsupported-shape dispatch:

- A polynomial matrix with determinant not equal to `1` throws an
  `ArgumentError` that mentions the determinant/unit precondition and says the
  input is outside the staged `SL_n` path.

Laurent determinant errors are delegated to `normalize_laurent_gl_matrix(A)`:

- Non-unit and unsupported unit determinants keep the staged messages from
  issue #20.

Valid but unsupported inputs reach algorithm-layer errors:

- Larger determinant-one polynomial inputs mention that the
  `SL_n` reduction layer for `n > 3` is not implemented.
- Multivariate or non-local `3 x 3` polynomial inputs mention the missing
  staged reduction to the local univariate `SL_3` slice.
- Normalized Laurent inputs mention that Laurent `SL_n` reduction is not
  implemented after the `GL_n` normalization boundary.

## Tests

Add `test/public/factorization_driver_shell.jl` and register it in the public
test group. The focused test verifies:

- The existing supported `3 x 3` univariate local `SL_3` example still
  factorizes exactly.
- A valid larger determinant-one polynomial matrix reaches a staged
  `SL_n` reduction-layer error, not a `3 x 3` hard stop.
- Nonsquare input is rejected by generic validation.
- A square polynomial determinant-not-one input fails on the determinant/unit
  precondition, not unsupported shape.
- A normalizable Laurent `GL_n` input reaches a post-normalization Laurent
  algorithm-layer error.
- A non-normalizable Laurent determinant keeps the issue #20 determinant
  precondition error.

The required focused verification command is:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Final verification uses:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl all
```

## Spec Self-Review

- No placeholders or incomplete markers remain.
- The design is scoped to the public factorization driver shell and focused
  public tests.
- The design preserves the existing supported factor sequence output.
- Determinant normalization remains delegated to the issue #20 boundary.
- The negative control cannot fall through to an unsupported-shape error
  because determinant checks happen before staged algorithm dispatch.
