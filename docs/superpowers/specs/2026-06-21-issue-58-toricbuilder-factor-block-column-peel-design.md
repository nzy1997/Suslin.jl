# Issue 58 ToricBuilder Factor Block Column-Peel Design

## Goal

Promote the checked-in ToricBuilder `factor_toric_block_3_pinv` and
`factor_toric_block_3_qinv` determinant-one fixtures from contract-only records
to exact elementary factorizations through the recursive Laurent column-peel
route.

## Context

Issue #57 merged the recursive Laurent column-peel reducer and routes Laurent
determinant-one inputs through it after the existing block-local `SL_n -> SL_3`
path. The ToricBuilder factor-block fixtures are larger determinant-one Laurent
matrices over `GF(2)[x^+/-1, y^+/-1]`:

- `factor_toric_block_3_pinv`: `8 x 8`, final peel block `I_2`.
- `factor_toric_block_3_qinv`: `16 x 16`, final peel block
  `[y 0; x*y^-1 + y^-1 y^-1]`.

A local probe on this branch confirms the existing reducer already reaches
those final blocks and produces exact nonempty factor sequences for both.
Therefore this issue should add durable public acceptance coverage and update
shared catalog metadata, not broaden the algorithm beyond this fixture family.

## Approaches Considered

1. Add focused public acceptance coverage and catalog promotion for Pinv/Qinv.

   This is the selected approach. It matches the issue objective, uses only the
   checked-in fixtures, and verifies the exact column-peel replay before
   changing shared catalog status.

2. Add fixture-specific branches inside the column-peel algorithm.

   This is unnecessary because the generic #57 determinant-one path already
   factors both fixtures exactly. A fixture-specific branch would add risk
   without adding mathematical coverage.

3. Depend on ToricBuilder at test time to regenerate the matrices.

   This is out of scope. The repository intentionally records offline fixtures
   and must not add a ToricBuilder runtime dependency.

## Design

Create `test/public/toricbuilder_factor_toric_block_acceptance.jl`. The test
will load `test/fixtures/toricbuilder_factor_toric_block_3.jl`, locate the
`Pinv` and `Qinv` entries, and assert:

- the fixture inverse relation remains exact;
- `elementary_factorization(entry.matrix)` returns a nonempty factor sequence;
- `verify_factorization(entry.matrix, factors)` is true;
- `Suslin._factor_laurent_sl_column_peel(entry.matrix)` verifies its replay;
- the replay dimensions are `[8, 7, 6, 5, 4, 3]` for `Pinv` and
  `[16, 15, ..., 3]` for `Qinv`;
- the final `2 x 2` blocks are exactly `I_2` for `Pinv` and
  `[y 0; x*y^-1 + y^-1 y^-1]` for `Qinv`;
- replacing one returned factor with identity fails exact verification;
- corrupting the stored ToricBuilder inverse relation fails the relation check.

Update `test/fixtures/toricbuilder_laurent_problem_catalog.jl` so the two
contract entries record `:supported_column_peel`, the `:laurent_column_peel`
path, the public acceptance verifier, and issue #58 as a consumer. Update
`test/internal/toricbuilder_problem_catalog.jl` so catalog validation accepts
and asserts that promoted status.

Register the new public test in `test/runtests.jl`.

## Error Handling

The public acceptance test will fail loudly through normal `@test` assertions
if the fixture relation is corrupted, if the factor sequence is empty, if exact
factor verification fails, or if the replay final block changes. The catalog
validator will reject unsupported status values, missing verifier paths, and
missing metadata exactly as it does for existing entries.

## Testing

Focused verification:

```bash
julia --project=. -e 'include("test/public/toricbuilder_factor_toric_block_acceptance.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_problem_catalog.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The acceptance file itself is the issue's required public command. It also
contains the negative controls requested by the issue.

## Scope Boundaries

- Do not add a ToricBuilder dependency.
- Do not require local ToricBuilder cache examples to pass.
- Do not support arbitrary determinant-one Laurent matrices beyond the
  existing column-peel behavior.
- Do not add or commit `Manifest.toml`.
