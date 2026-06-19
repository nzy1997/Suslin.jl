# Issue 17 Laurent Large Acceptance Design

## Context

Issue #17 is the end-to-end acceptance gate for the Laurent factorization
roadmap. The dependency issues have already added the ToricBuilder contract
fixtures, the Laurent determinant-normalization boundary, the public driver
shell, local `SL_3` solving, and the first `SL_n` to local `SL_3` reduction
record. This issue should not introduce an unrelated general Suslin
algorithm. It should verify the landed supported path at the target matrix
scale and keep unsupported real-consumer inputs visible.

The current ToricBuilder contract from issue #19 records real
`factor_toric_block(3, x, y, R)` matrices over
`GF(2)[x^+/-1, y^+/-1]`. Those fixtures remain larger two-variable Laurent
inputs whose current contract is a verified transformation certificate and
normalization boundary, not a complete elementary factorization.

## Approaches Considered

1. Add a dedicated acceptance catalog plus a public acceptance test. This is
   the selected approach because it keeps the large fixtures reusable, records
   provenance near the matrix constructors, and verifies the exact public
   driver behavior in one focused gate.
2. Put all generated matrices inline inside `test/public`. That would satisfy
   the test mechanically, but it would make the 40x40 and larger examples hard
   to inspect or reuse.
3. Expand to arbitrary Laurent matrix factorization. That is broader than an
   acceptance harness and would blur staged unsupported cases that the
   dependency issues deliberately preserved.

## Selected Design

Create `test/fixtures/laurent_large_acceptance_cases.jl` with a small catalog
module. The module will expose `acceptance_catalog()` returning cases with:

- `id`
- `kind`
- `ring` metadata
- `matrix`
- `size`
- `provenance`
- `expected_path`
- `negative_control`

The catalog contains three acceptance rungs:

1. A real ToricBuilder-motivated fixture from issue #19, using the checked-in
   `Pinv` relation from `factor_toric_block(3, x, y, R)`. This case follows
   the normalized/contract path: `normalize_laurent_gl_matrix(A)` and
   `verify_laurent_gl_normalization(A, normalization)` must succeed, and the
   stored ToricBuilder inverse relation must still multiply exactly to the
   identity.
2. A 40x40 Laurent `SL_n` matrix over `GF(2)[x^+/-1, y^+/-1]` built from
   disjoint embedded local `SL_3` blocks. This case must factorize through
   `elementary_factorization(A)` and verify with `verify_factorization(A,
   factors)`.
3. A 48x48 Laurent `SL_n` matrix with the same supported block-local shape.
   Runtime probing in this Agent Desk worktree showed determinant
   normalization for the generated 48x48 case completes in a few seconds, so
   this larger case is enabled by default.

The large matrices are synthetic, but they do not replace the real consumer
fixture. They exist only to exercise the supported exact factorization path at
the required scale.

## Minimal Implementation Needed

The landed `SL_n` reduction design already describes Laurent normalized
cores: Laurent inputs should pass through determinant normalization, use the
first Laurent generator for local obligations, and call `realize_sl3_local`
with `check_monic=false`. The current implementation still stages every
Laurent reduction as unsupported. Issue #17 should finish that narrow supported
path only for block-local Laurent determinant-one cores:

- permit Laurent rings in `reduce_sln_to_sl3` instead of throwing in
  `_reduction_generator`,
- pass `check_monic=false` for Laurent local obligations,
- record Laurent-specific assumptions on each obligation,
- keep `elementary_factorization(A)` returning factors only when the assembled
  factors exactly multiply to the original input.

If a Laurent `GL_n` case has a non-identity determinant correction, the driver
still must not claim a raw elementary factorization of the original matrix
unless exact factor reconstruction verifies.

## Acceptance Test

Create `test/public/laurent_large_acceptance.jl`. The test will:

- include the large acceptance catalog module,
- assert that the catalog contains one ToricBuilder case, one 40x40 case, and
  one enabled larger case,
- run the normalized ToricBuilder path and exact inverse-relation verifier,
- run `elementary_factorization(A)` for every factorization case,
- assert `verify_factorization(A, factors)` and `verify_sln_to_sl3_reduction`
  for the reduction metadata,
- include a negative control that corrupts one returned factor and proves exact
  verification fails.

Register the new public test in `test/runtests.jl` so the package entry point
and the documented full suite both cover the acceptance gate.

## Verification

Run the issue-specific command:

```bash
julia --project=. -e 'include("test/public/laurent_large_acceptance.jl")'
```

Run the Agent Desk package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run the documented full suite from issue #21:

```bash
julia --project=. test/runtests.jl all
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design and written spec are approved automatically because the selected
approach is the narrowest exact acceptance gate that satisfies issue #17,
keeps the real ToricBuilder fixture in scope, and limits algorithm changes to
the already-designed Laurent block-local reduction path.

## Spec Self-Review

- No incomplete markers remain.
- The ToricBuilder fixture cannot be replaced by synthetic large matrices.
- The 40x40 and 48x48 cases exercise exact factorization and exact
  verification.
- The negative control fails by exact verification.
- The scope excludes arbitrary Laurent factorization and non-elementary
  determinant-correction claims.
