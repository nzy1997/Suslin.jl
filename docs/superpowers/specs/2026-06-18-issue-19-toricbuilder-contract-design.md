# Issue 19 ToricBuilder Contract Design

## Context

Suslin currently factors only a narrow univariate polynomial `SL_3` slice.
ToricBuilder's consumer boundary is different: it produces exact GF(2)
two-variable Laurent matrices such as `Qinv`, `Pinv`, `row_transformation`,
and `column_transformation`. The first integration contract should record
small real ToricBuilder cases and their expected Suslin status before broader
factorization work expands the accepted input class.

The local ToricBuilder checkout used for provenance is:

- Path: `/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder`
- Commit: `fa7f82252d42fdc0b2726bc48af24ac4c70a8d73`
- Source function: `src/toric_form/toric_factorization.jl:factor_toric_block`
- Generation entry point: `factor_toric_block(3, x, y, R)` with
  `R, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])`

## Design Choice

Use a checked-in Julia fixture plus a docs contract.

The fixture will live under `test/fixtures/` and contain literal Oscar Laurent
matrix expressions and metadata. Tests will include the fixture and verify the
metadata, determinant classification, and exact multiplication relations.
Runtime Suslin code will not depend on ToricBuilder.

Alternative considered: add only an optional regeneration script. That keeps
fixtures fresher, but it makes the core contract hard to test on machines that
do not have ToricBuilder available.

Alternative considered: add a public Suslin fixture-loading API. That is too
large for this issue because the fixture is a test contract, not a supported
library interface.

## Contract Scope

The first supported contract entries are the two inverse-style matrices from
`factor_toric_block(3, x, y, R)`:

- `Qinv`: a `16 x 16` matrix over `GF(2)[x^+/-1, y^+/-1]`.
- `Pinv`: an `8 x 8` matrix over `GF(2)[x^+/-1, y^+/-1]`.

Both entries record determinant classification `one`. The fixture also stores
their ToricBuilder source matrices:

- `column_transformation * Qinv == I_16`
- `row_transformation * Pinv == I_8`

The fixture records the block coarse-graining size `(2, 2)` returned by
`factor_toric_block(3, x, y, R)`.

## Expected Suslin Behavior

Suslin is not expected to factor these matrices yet. Current behavior should
remain staged unsupported input:

- `elementary_factorization(Qinv)` fails before factorization because the
  matrix is not `3 x 3` and is over a two-variable Laurent ring.
- `elementary_factorization(Pinv)` fails before factorization for the same
  class of reasons.

Later implementation issues can decide whether to normalize Laurent monomial
determinants before factorization. This issue only records the real boundary.

## ToricBuilder Output Need

The downstream need is a verified transformation certificate, not just raw
elementary factors. For this contract, the certificate is the exact inverse
relation above plus determinant classification. Later factorization output can
add raw elementary factors and normalized metadata, but those are not required
for this issue.

## Files

- Create `docs/src/toricbuilder_contract.md`: human-readable integration
  contract, fixture provenance, determinant classification rules, and expected
  current behavior.
- Modify `docs/src/index.md`: link the new contract page.
- Create `test/fixtures/toricbuilder_factor_toric_block_3.jl`: literal fixture
  data and metadata for `Qinv` and `Pinv`.
- Create `test/internal/toricbuilder_contract.jl`: exact fixture checks and
  negative controls.
- Modify `test/runtests.jl`: include the new internal test file.

## Verification

The focused test must verify:

- At least one `Qinv` and one `Pinv`-style matrix are represented exactly.
- Every fixture records ring, size, determinant classification, expected
  behavior, and provenance.
- Each claimed ToricBuilder inverse relation is checked by exact multiplication.
- A deliberately wrong determinant classification fails fixture validation.
- A deliberately corrupted matrix relation fails fixture validation.

The final verification command is the issue 21 package entry point:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The full suite command from issue 21 is also useful as a broader check:

```bash
julia --project=. test/runtests.jl all
```

## Spec Self-Review

- No incomplete markers remain.
- The scope is limited to contract docs, checked-in fixture data, and tests.
- The design avoids public API changes and full factorization work.
- The documented expected behavior matches Suslin's current implementation
  limits.
