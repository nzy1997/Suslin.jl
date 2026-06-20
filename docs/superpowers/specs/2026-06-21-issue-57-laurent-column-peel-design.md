# Issue 57 Laurent Column-Peel Design

## Goal

Add a deterministic recursive Laurent `SL_n` column-peel factorization path for
the Issue #38 row-normalized and column-normalized determinant-one cores, and
route those cores through `elementary_factorization`.

## Context

Issue #38 records a `6 x 6` Laurent `Q` block over
`GF(2)[u^+/-1, v^+/-1]` with determinant `u*v`. The original `Q` remains a
Laurent `GL_n` input and still needs a later determinant-correction
certificate. Issue #57 is narrower: after row or column Laurent-unit
normalization, the determinant-one cores should factor into elementary
matrices.

The existing `reduce_sln_to_sl3` path only handles disjoint embedded local
`SL_3` obligations. Diagnostics from Issue #41 show the Issue #38 row and
column cores fail that block-local shape. A direct probe on this branch showed
the last column is reducible at each stage by `reduce_unimodular_column`; after
left column reduction and right bottom-row clearing, the path reaches these
final blocks exactly:

```text
row core final 2x2:    [u  u*v; 0  u^-1]
column core final 2x2: [v^-1  u*v; 0  v]
```

## Approaches Considered

1. Add a focused Laurent column-peel reducer and use it as the Laurent
   determinant-one fallback in `elementary_factorization`.

   This is the selected approach. It follows the issue's Park-Woodburn-aligned
   algorithm, reuses existing exact column reduction and local `SL_3` solver
   machinery, and keeps unsupported cases explicit.

2. Extend the existing `SL_n -> local SL_3` block-location search.

   This would keep one reduction path, but Issue #41 already showed that no
   `3+3` partition solves these cores. Extending that search is the wrong
   bottleneck for this issue.

3. Add a broad matrix search for elementary factors.

   This is out of scope. It would be harder to make deterministic and would not
   exercise the specific recursive column-peel bridge requested by Issue #57.

## Design

Create `src/algorithm/laurent_column_peel.jl` and include it from
`src/Suslin.jl` after the local `SL_3` solver is available. The implementation
adds internal replay metadata types:

- `LaurentColumnPeelStep`: one peel level, including dimension, input matrix,
  last column, left factors, after-left matrix, right bottom-row clearing
  factors, peeled matrix, and upper-left next block.
- `LaurentColumnPeelFactorization`: original matrix, final `2 x 2` block, final
  local factors, full elementary factors, product, peel steps, and verification
  data.

The main internal entry point is `_factor_laurent_sl_column_peel(A)`. It accepts
only square Laurent determinant-one matrices of size at least `2`. Public
`elementary_factorization(A)` will keep trying the existing block-local
`reduce_sln_to_sl3` path first for Laurent determinant-one inputs, then fall
back to `_factor_laurent_sl_column_peel(A)` when the block-local path throws a
staged `ArgumentError`.

Each recursive peel step does this on the current `d x d` matrix:

1. Take column `d`.
2. Call `reduce_unimodular_column` on that column. Its factors multiply on the
   left and send the column to `e_d`.
3. Build right elementary matrices `E(d, j, -B[d, j])` for `j < d` after the
   left reduction. Multiplying by them clears the bottom row because the last
   column is `e_d`.
4. Verify the peeled matrix is block diagonal with upper-left `(d-1) x (d-1)`
   block and bottom-right `1`.
5. Recurse on the upper-left block.

If `L * current * C = blockdiag(next, 1)`, and the recursive factorization of
`next` is embedded as `P`, the current matrix factors as:

```text
current = inv(L) * P * inv(C)
```

The returned factor sequence therefore appends the inverse left factors in
reverse order, the embedded recursive factors, and the inverse right-clearing
factors in reverse order.

The final `2 x 2` determinant-one block is embedded as a `3 x 3` local matrix
with trailing identity and solved with `realize_sl3_local(...; check_monic =
false)`. Because the supported unit-pivot local solver emits only `E(1,2,*)`
and `E(2,1,*)` factors for these final Issue #38 blocks, the trailing coordinate
is dropped back to `2 x 2` elementary factors.

## Verification

Add `test/expert/laurent_column_peel_issue38.jl` and register it in the expert
group. The test should:

- Load `test/fixtures/toricbuilder_issue38_cases.jl`.
- Factorize row and column determinant-one cores through the column-peel helper
  and through `elementary_factorization`.
- Assert the factor sequences are nonempty and
  `verify_factorization(core, factors)` is true.
- Assert replay metadata dimensions are `6, 5, 4, 3`.
- Assert every recorded peel step satisfies exact
  `left_product * input * right_product == peeled_matrix` reconstruction and
  exposes the next upper-left block.
- Assert the final `2 x 2` blocks equal
  `[u  u*v; 0  u^-1]` for the row core and `[v^-1  u*v; 0  v]` for the column
  core.
- Corrupt one recorded peel left factor and delete one bottom-row clearing
  factor; both replay verification and final `verify_factorization` checks must
  fail.

Update `test/internal/toricbuilder_issue38_fixture.jl` so the fixture still
validates the original `Q` determinant metadata and normalization data, but now
records row and column normalized cores as supported by the column-peel path
instead of expected `elementary_factorization` failures. Keep the original `Q`
factorization boundary unsupported.

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_peel_issue38.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not factorize the original Issue #38 `Q` matrix with determinant `u*v`. Do
not add a general determinant-correction certificate. Do not add broad matrix
search or support arbitrary Laurent `SL_n` inputs beyond the exact staged
column-peel cases accepted by existing column reduction and final unit-pivot
local `SL_3`.

## Spec Self-Review

- No placeholders or incomplete sections remain.
- The design implements the Issue #57 algorithm directly and keeps the original
  `GL_n` determinant-correction problem out of scope.
- The fallback order preserves currently supported block-local Laurent cases
  while adding the Issue #38 determinant-one core path.
- The replay metadata has enough exact matrices and factors to verify each peel
  step and negative controls.
