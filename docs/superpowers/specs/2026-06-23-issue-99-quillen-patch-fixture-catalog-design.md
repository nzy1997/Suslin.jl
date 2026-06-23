# Issue 99 Quillen Patch Fixture Catalog Design

## Context

Issue 63 will need deterministic local-to-global Quillen patching examples
before the patching layer grows more theorem-bearing API. The repository
already has ordinary polynomial Quillen examples in
`test/expert/quillen_patching_exact.jl`, plus a separate patched-substitution
scaffold in `test/expert/quillen_induction.jl`. Later replay issues need these
examples as named fixtures with exact denominator coverage, local elementary
corrections, expected product data, and explicit patched-substitution witness
fields.

Issue 99 asks for test support only. Laurent `GL_n`, ToricBuilder-derived
examples, constructive patch assembly, and public API changes stay out of
scope.

## Design Choice

Add a test-only fixture module at `test/fixtures/quillen_patch_cases.jl` and a
focused internal validator at `test/internal/quillen_patch_fixtures.jl`.

Each fixture is a named tuple with:

- `id`
- `kind`
- `stage_coverage`
- `ring_constructor`
- `ring`
- `size`
- `substitution_variable`
- `target_matrix`
- `base_matrix`
- `denominator_data`
- `local_factors`
- `expected`
- `patched_substitution_witness`
- `source_refs`
- `consumer_issue_ids`

The catalog returns valid cases and negative controls. Negative controls are
deliberately corrupted entries that must fail the same validator as valid
entries, proving coverage multipliers and local factors are checked as data
with algebraic meaning rather than decorative metadata.

Alternatives considered:

- Promote fixtures through exported Suslin APIs. That would add compatibility
  burden before Issue 63 chooses a constructive API.
- Keep examples only inside `test/expert/quillen_patching_exact.jl`. That
  would keep tests green but would not give later issues a shared catalog.
- Build a generic random example bank. The issue asks for deterministic exact
  ordinary-polynomial cases, so the catalog stays small and explicit.

## Fixture Scope

The valid catalog contains at least five entries:

- `quillen-two-open-cover-qq`: seeded from the ordinary `QQ[X,r,g]` two-open
  cover in the current expert test.
- `quillen-nontrivial-multipliers-qq`: an ordinary `QQ[X,r,g]` two-open cover
  whose Bezout coverage multipliers are nonconstant rather than both `1`.
- `quillen-supplied-local-certificate-gf2`: an ordinary `GF(2)[X,r,s]`
  example with supplied local certificates and exact local factors.
- `quillen-patched-substitution-witness-qq`: a `QQ[X,r,g]` example that records
  `X`, denominator `r`, exponent `l`, and shift `g` for later
  `patched_substitution(A, X, r, l, g)` replay.
- `quillen-constructive-acceptance-gf2`: a second constructive ordinary
  polynomial acceptance case over `GF(2)[X,r,s]`.

The negative controls include:

- `quillen-uncovered-denominator-control`: mutates one coverage multiplier so
  `sum(c_i * d_i) != 1`.
- `quillen-tampered-local-factor-control`: mutates one local factor so the
  local factor product no longer equals the expected correction.

No Laurent examples are migrated into this catalog.

## Validation Rules

The validator checks required metadata, unique ids, ordinary polynomial ring
metadata, matrix sizes, generator-valued substitution variables, denominator
coverage equations, local certificate shape, and local factor products.

For positive entries it reconstructs:

- the ring and named generators,
- target and optional base matrices,
- denominator and coverage multiplier pairs,
- supplied local certificates,
- local elementary factor matrices,
- expected global correction,
- patched-substitution witness values where present.

It checks `sum(c_i * d_i) == 1` exactly, each local factor equals the elementary
matrix described by its local correction, the product of local factors equals
the expected correction, and determinant-one claims hold exactly. Fixtures with
`expected.current_status == :passes` must validate; fixtures with
`expected.current_status == :staged_fail` are allowed only when they carry
explicit missing witness metadata.

Patched-substitution witnesses must include `matrix`, `variable`, `denominator`,
`exponent`, `shift`, and `expected_matrix`. The validator calls
`Suslin.patched_substitution(matrix, variable, denominator, exponent, shift)`
and compares the result to `expected_matrix`, so Issue 102 can later replay and
tamper-check these fields without changing the schema.

The catalog-level validator also validates negative controls by requiring them
to throw `ArgumentError`. The internal test mutates one valid coverage
multiplier and one valid local factor at runtime and verifies both mutations are
rejected.

## Files

- Create `test/fixtures/quillen_patch_cases.jl`: exact ordinary polynomial
  fixture construction and helper lookup by id.
- Create `test/internal/quillen_patch_fixtures.jl`: validator functions,
  focused tests, id coverage checks, and negative controls.
- Modify `test/runtests.jl`: include the validator in the `internal` group.

## Verification

Focused validator:

```bash
julia --project=. -e 'include("test/internal/quillen_patch_fixtures.jl")'
```

Internal group:

```bash
julia --project=. test/runtests.jl internal
```

Full package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to test fixtures and internal validation.
- Required positive and negative fixture ids are listed.
- Laurent and ToricBuilder cases remain out of scope.
- No public Suslin API changes are proposed.
