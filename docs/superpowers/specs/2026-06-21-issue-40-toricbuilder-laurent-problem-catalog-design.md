# Issue 40 ToricBuilder Laurent Problem Catalog Design

## Context

Issue #40 asks for a single reusable manifest of ToricBuilder-related Laurent
problem cases. The relevant inputs already exist as test fixtures:

- `test/fixtures/toricbuilder_issue38_cases.jl` contains the Issue #38 `6 x 6`
  `Q` block and its row and column determinant-one cores.
- `test/fixtures/toricbuilder_factor_toric_block_3.jl` contains the checked-in
  `factor_toric_block_3_qinv` and `factor_toric_block_3_pinv` contract
  matrices from the ToricBuilder boundary.
- `test/fixtures/laurent_large_acceptance_cases.jl` contains the synthetic
  `40 x 40` and `48 x 48` block-local Laurent acceptance cases.

Later algorithm tests should consume one catalog entry per case instead of
redefining stable IDs, provenance, determinant metadata, expected current
status, verifier paths, or milestone links.

## Design Choice

Add a dedicated manifest-style catalog at
`test/fixtures/toricbuilder_laurent_problem_catalog.jl` and a focused internal
validator at `test/internal/toricbuilder_problem_catalog.jl`.

This is the safest shape because the catalog spans multiple existing fixture
modules and has different metadata than the older general Laurent fixture
catalog. A dedicated ToricBuilder problem catalog avoids changing public APIs,
does not require ToricBuilder at runtime, and keeps consuming tests pointed at a
single source of truth.

Alternatives considered:

- Extend `test/fixtures/laurent_cases.jl`. That catalog is organized around
  low-level Laurent algebra relations and would need unrelated fields for
  expected ToricBuilder algorithm status and verifier routing.
- Put the manifest directly in the internal test. That would validate the
  entries but would not provide a reusable fixture file for public, internal,
  and expert tests.

## Catalog Interface

The fixture module exports a module-level function:

```julia
ToricBuilderLaurentProblemCatalog.catalog()
```

It returns a named tuple with `cases`, where each entry has these required
fields:

- `id`: stable string ID.
- `kind`: symbolic family such as `:issue38_q_block`,
  `:toricbuilder_contract`, or `:synthetic_block_local_acceptance`.
- `source_fixture`: where the matrix object comes from.
- `ring`: description plus object/generator metadata when available.
- `dimensions`: matrix size metadata.
- `matrix`: the matrix consumed by verifiers.
- `determinant_profile`: determinant class and expected determinant metadata.
- `expected_current_status`: status symbol such as `:unsupported_now`,
  `:verified_contract`, `:supported_block_local`, or `:target_acceptance`.
- `verifier`: verifier path and scenario metadata.
- `provenance`: issue/source data from the originating fixture.
- `consumers`: milestone and issue metadata for current or planned consuming
  tests.

The initial catalog contains five entries:

- `toricbuilder-issue-38-q-block` with determinant class
  `:laurent_monomial_unit`, status `:unsupported_now`, and verifier path
  `test/internal/toricbuilder_issue38_fixture.jl`.
- `toricbuilder-factor-toric-block-3-qinv` with determinant class `:one`,
  status `:verified_contract`, and verifier path
  `test/internal/toricbuilder_contract.jl`.
- `toricbuilder-factor-toric-block-3-pinv` with determinant class `:one`,
  status `:verified_contract`, and verifier path
  `test/internal/toricbuilder_contract.jl`.
- `laurent-block-local-40x40` with determinant class `:one`, status
  `:supported_block_local`, and verifier path
  `test/public/laurent_large_acceptance.jl`.
- `laurent-block-local-48x48` with determinant class `:one`, status
  `:target_acceptance`, and verifier path
  `test/public/laurent_large_acceptance.jl`.

The `48 x 48` synthetic case already passes through the same acceptance test as
the `40 x 40` case. The `:target_acceptance` status records that it is the
larger future-facing acceptance representative named by Issue #40, while
retaining its concrete verifier path.

## Validation

The internal validator checks the catalog as metadata and as live fixture
references:

- IDs are unique.
- Every entry has nonempty provenance, determinant metadata, expected status,
  verifier path, and consumer metadata.
- Required IDs from Issue #40 are present.
- The recorded matrix dimensions match the actual matrix.
- Determinant metadata matches the actual determinant classification.
- Verifier paths point to checked-in test files.

Negative controls construct catalog variants with a duplicate ID and missing
provenance. Both must be rejected by the same validator, proving the manifest is
not accepted by shape alone.

## Out of Scope

This issue does not implement new factorization algorithms, change existing
ToricBuilder fixtures, or require a local ToricBuilder checkout. It only adds a
shared manifest and validation tests.

## Verification

Focused verification:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_problem_catalog.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No placeholders or incomplete sections remain.
- The design is scoped to one fixture module, one internal validator, and test
  registration.
- The catalog consumes existing fixtures instead of duplicating matrices.
- The negative controls have exact observable failure modes.
