# Issue 154 Lazy Laurent Determinant Fixture Catalog Design

## Context

Issue #154 prepares later lazy Laurent determinant work. The current Laurent
`GL_n` route normalizes an input matrix by computing its full determinant up
front. Large ToricBuilder Q-block runs can spend too much time in that
determinant and normalization work. This issue does not implement lazy
determinant peeling; it adds a shared offline catalog of compact fixtures that
future lazy-determinant issues can cite.

The repository already has a source-grounded Issue #38 ToricBuilder fixture at
`test/fixtures/toricbuilder_issue38_cases.jl`. That fixture constructs the
`6 x 6` Q block over `GF(2)[u^+/-1, v^+/-1]`, records determinant `u*v`, and
stores row and column determinant-one cores. The new catalog should wrap that
entry instead of copying the matrix.

## Approach Options

Recommended: add a dedicated fixture module at
`test/fixtures/laurent_lazy_determinant_cases.jl` and a dedicated internal
validator at `test/internal/laurent_lazy_determinant_fixtures.jl`. The module
returns named cases with a consistent lazy-determinant metadata shape. The
validator performs exact determinant, correction-support, provenance, and
negative-control checks. This follows the existing fixture-catalog pattern and
keeps the new work out of public APIs.

Alternative: append the cases to `test/fixtures/laurent_cases.jl`. This would
reuse the general Laurent catalog, but the lazy determinant cases need row and
column normalized-core metadata plus correction-support expectations that do
not fit the older catalog shape.

Alternative: add only an internal test without a fixture module. This would
exercise the expected behavior, but it would not give later issues a stable
catalog to cite.

## Chosen Design

Add module `LaurentLazyDeterminantCases` with `catalog()` returning four cases:

- `issue-38-q-block-lazy-determinant`: wraps
  `ToricBuilderIssue38Cases.catalog()` and records determinant `u*v`,
  classification `:laurent_monomial_unit`, supported correction metadata, and
  the existing row and column determinant-one cores.
- `determinant-one-triangular`: a compact synthetic determinant-one Laurent
  matrix that requires no correction.
- `monomial-unit-row-column-cores`: a compact synthetic Laurent monomial-unit
  matrix whose row and column normalized cores both have determinant one.
- `non-unit-determinant-negative`: a compact synthetic matrix with determinant
  `x + 1`, expected classification `:non_unit`, and unsupported correction
  metadata.

Each entry records:

- `id` and `kind`;
- `ring` metadata with the Oscar parent, generators, variables, and
  description;
- `dimensions.matrix`;
- `inputs.matrix`;
- `determinant_profile.expected_determinant`,
  `determinant_profile.expected_class`, and monomial metadata where relevant;
- `expected_correction.supported`, `expected_correction.kind`, and
  `expected_correction.supports`;
- `normalizations.row` and `normalizations.column` for supported cases;
- `negative_control` metadata describing how the validator should reject drift;
- `provenance` with source and issue fields;
- `consumer_test_ids`.

For supported monomial-unit cases, the row core comes from
`normalize_laurent_gl_matrix(A)` and the column core uses a right diagonal
factor `diag(inv(det(A)), 1, ...)`. For determinant-one cases, the row and
column cores are the original matrix and both corrections are identities. For
the non-unit case, normalizations are omitted and the validator asserts that
the supported-correction predicate rejects it.

## Validation

The internal validator exposes:

```julia
validate_laurent_lazy_determinant_fixture(entry)
validate_laurent_lazy_determinant_catalog(catalog)
_fixture_supports_lazy_determinant_correction(entry)
```

Validation checks required fields, unique IDs, exact dimensions, determinant
metadata, correction support, provenance fields, and consumer-test IDs. It
recomputes `classify_laurent_determinant(entry.inputs.matrix)` and rejects any
entry whose expected determinant class or determinant value drifts from the
matrix.

For supported entries, the validator checks the row and column normalized
cores exactly:

- row normalization verifies with `verify_laurent_gl_normalization`;
- column normalization reconstructs as `A * correction_factor`;
- both cores have determinant one;
- where practical, determinant-one cores are also passed through existing
  factorization helpers or exact determinant checks.

The Issue #38 wrapper has an extra drift check against
`ToricBuilderIssue38Cases.catalog()`: the wrapped matrix, determinant `u*v`,
row core, and column core must match the source fixture exactly.

The validator includes two required negative controls:

- change a monomial-unit fixture's expected class to `:one` and assert that
  validation throws `ArgumentError`;
- assert the non-unit fixture is rejected by
  `_fixture_supports_lazy_determinant_correction(entry)` instead of being
  marked correctable.

## Test Registration

Register `internal/laurent_lazy_determinant_fixtures.jl` in the `internal`
test group in `test/runtests.jl`.

Focused verification:

```bash
julia --project=. -e 'include("test/internal/laurent_lazy_determinant_fixtures.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not implement lazy determinant peeling. Do not change public APIs. Do not
change the current Laurent `GL_n` normalization behavior. Do not duplicate the
Issue #38 matrix outside the existing Issue #38 fixture module.

## Automatic Decisions

- Visual companion: skipped because this is a metadata and algebraic validation
  task, not a visual design task.
- Clarifying questions: resolved from the issue text because Agent Desk marked
  this run non-interactive.
- Approach: choose the dedicated catalog and validator because it matches
  existing fixture patterns and keeps the behavioral surface internal.
- Issue #38 source: wrap `ToricBuilderIssue38Cases.catalog()` to avoid matrix
  drift and satisfy the issue's reuse requirement.

## Spec Self-Review

- No placeholders or incomplete sections remain.
- The scope is limited to a fixture module, validator test, runner
  registration, and workflow docs.
- The fixture fields cover the issue's required matrix, ring, determinant,
  correction support, and negative-control metadata.
- The validation strategy has exact observable failures for both catalog drift
  and unsupported non-unit correction.
