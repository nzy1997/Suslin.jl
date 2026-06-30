# Issue 217 Park-Woodburn Substitution Chain Design

## Context

Issue #217 follows the #216 bounded denominator-cover solver. The repository
now has `QuillenDenominatorCoverSolverResult`, which records raw denominators
`r_i`, a chosen exponent `l`, powered terms `r_i^l`, multipliers `g_i`, exact
coverage terms, and a replayable cover certificate proving
`sum(g_i * r_i^l) == 1`.

Park-Woodburn Section 3 needs the next internal certificate layer: a replayable
chain of cumulative substitutions
`X -> X - X * sum_{j <= i} g_j * r_j^l`. The final cumulative coefficient is
zero when the cover identity holds, the final base term is `A(0)`, and the
bracket terms telescope from `A(X)` to `A(0)`. This chain is separate from the
older `patched_substitution(A, X, r, l, g)` helper, which records the older
`X + r^l*g` witness shape.

## Approaches Considered

Recommended: add an internal Park-Woodburn chain certificate in
`src/algorithm/quillen_induction.jl`, backed by a small matrix substitution
helper that evaluates `A` at `X => coefficient * X`. Construction consumes a
verified #216 solver result, records each coefficient and matrix, records each
bracket `A(previous)^{-1} * A(next)`, and verifies the telescope by exact
matrix multiplication. This keeps all Quillen replay records together and
reuses the existing exact polynomial-ring validation patterns.

Alternative: extend the old `patched_substitution` helper to accept a
coefficient mode. This would blur two different witness shapes and make later
global factor assembly infer Section 3 metadata from an older interface.

Alternative: store only coefficient metadata and ask consumers to recompute
matrices. That is smaller, but it would not meet the issue objective because
the certificate itself would not reject tampered intermediates, brackets, or
base terms.

The chosen design is the recommended internal chain certificate because it is
replayable, exact, narrowly scoped, and leaves public factorization behavior
unchanged.

## Certificate Shape

Add records beside the #216 solver records:

- `QuillenPatchSubstitutionStep`: one denominator-cover step. It records the
  step index, selected variable, raw denominator, exponent, powered
  denominator, multiplier, sign convention, previous and next cumulative
  coefficients, previous and next substituted matrices, bracket target, and
  replay metadata.
- `QuillenPatchSubstitutionChainVerification`: replay data for the whole
  chain, including solver validity, ring and variable checks, coefficient and
  matrix replay, final cover identity status, base-term status, bracket
  telescope product, and overall status.
- `QuillenPatchSubstitutionChain`: the constructed certificate. It records the
  original matrix, ring, size, selected variable, sign convention, solver
  result, starting coefficient, final coefficient, intermediate matrices,
  steps, bracket matrices, base term `A(0)`, replay metadata, and verification.

The sign convention is explicit and fixed to `:park_woodburn_minus`, meaning
each step uses `next_coefficient = previous_coefficient - g_i * r_i^l`.
Because each substituted matrix is `A(coefficient * X)`, this represents
`X -> X - X * sum_{j <= i} g_j * r_j^l`.

## Construction And Replay

Add these internal functions:

- `_quillen_substitute_matrix_scaled_variable(A, X, coefficient)`: validates
  `X` as a generator of `base_ring(A)`, coerces the coefficient into the same
  ring, and evaluates each entry at `X => coefficient * X`.
- `quillen_patch_substitution_chain(A, X, solver_result; sign_convention = :park_woodburn_minus, metadata = (;))`.
- `replay_quillen_patch_substitution_chain(chain)`.
- `verify_quillen_patch_substitution_chain(chain)::Bool`.

Construction requires a square ordinary-polynomial matrix over the exact ring
stored in the solver result, a selected generator from that ring, and a solver
result passing `verify_quillen_denominator_cover_solver_result`. It starts from
coefficient `1`, then applies each powered denominator and multiplier in solver
order. It records `A(c_i X)` for every cumulative coefficient. The final
coefficient must be zero, so the final substituted matrix is `A(0)`.

Each bracket target is computed as
`inv(previous_matrix) * next_matrix`, using Oscar's exact matrix inverse. This
keeps the bracket direction aligned with the issue statement:
`A^{-1}(previous) * A(next)`. The telescope verifier multiplies
`A(X) * bracket_1 * ... * bracket_n` and checks that it equals the recorded
base term. This proves the recorded brackets telescope from the original
matrix to `A(0)`.

Replay recomputes every coefficient, substituted matrix, bracket, base term,
and telescope product from stored inputs and solver data. It compares the
stored chain and stored verification against that replay. Tampering with an
exponent, multiplier, sign convention, intermediate matrix, bracket, base
term, or selected variable makes verification return `false`.

## Tests

Create `test/expert/quillen_patch_substitution_chain.jl`. The positive case
uses the #213/#99 `quillen-two-open-cover-qq` fixture, solves its #216 cover,
builds a chain for `entry.base_matrix`, and proves:

- the result is a `QuillenPatchSubstitutionChain`;
- the solver and chain both verify;
- cumulative coefficients are `[1, 1 - r, 0]` for the two-open cover;
- intermediate matrices are exact evaluations of `A(cX)`;
- the base term equals exact `A(0)`;
- `A(X) * prod(brackets) == A(0)`;
- each step records selected variable, exponent, raw denominator, powered
  denominator, multiplier, sign convention, and replay metadata.

Negative controls rebuild the record with one corrupted exponent, multiplier,
sign convention, intermediate matrix, and selected variable. Each corrupted
record must be rejected by `verify_quillen_patch_substitution_chain`.

Register the focused expert test in `test/runtests.jl` after
`expert/quillen_denominator_cover_solver.jl`.

## Files

- Modify `src/algorithm/quillen_induction.jl`.
- Create `test/expert/quillen_patch_substitution_chain.jl`.
- Modify `test/runtests.jl`.
- Add `docs/superpowers/plans/2026-06-30-issue-217-park-woodburn-substitution-chain.md`.

Do not export the new names from `src/Suslin.jl`. Existing expert tests can
access unexported names through the `Suslin.` module qualifier. Do not factor
`A(0)`, solve local realizability, normalize local factors, assemble global
factors, or claim public `elementary_factorization` support.

## Verification

Focused chain verifier:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_substitution_chain.jl")'
```

Full package verifier:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Diff hygiene:

```bash
git diff --check
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it follows the issue's recommended
internal certificate shape, consumes the verified #216 solver result, records
the exact Park-Woodburn minus sign convention, and keeps out-of-scope
factorization work untouched.

## Spec Self-Review

- No placeholders remain.
- The design records every required coefficient, matrix, bracket target,
  selected variable, exponent, multiplier, sign convention, and replay field.
- The final coefficient and base term are tied directly to the cover identity.
- The telescope direction is explicit and tested.
- Negative controls cover all corruptions required by the issue.
- Scope is limited to the internal substitution-chain certificate.
