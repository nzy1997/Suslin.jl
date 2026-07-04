# Issue 288 Steinberg Optimization Fixture Catalog Design

Issue #288 adds a validated fixture catalog for Park-Woodburn Section 6
Steinberg rewrite opportunities. The catalog is test-only evidence for later
optimizer work and does not implement optimization or change public
factorization behavior.

## Context

There is no repository `AGENTS.md`; the README and existing fixture-catalog
patterns apply. Live `gh issue view 288` failed in this Agent Desk sandbox
because the configured GitHub proxy is not reachable, so the supplied issue
body, local branch state, Park-Woodburn Section 6 text in
`refs/arXiv-alg-geom9405003v1`, and existing test catalogs are the source of
truth. `gh pr list --search "288"` returned no matching PR context.

Park-Woodburn Section 6 lists the Steinberg relations:

- identity elementary factors `E_ij(0) = I`;
- same-position merges `E_ij(a)E_ij(b) = E_ij(a + b)`;
- two nontrivial commutator rewrites;
- the disjoint commutator identity.

The requested output is a shared fixture catalog for ordinary-polynomial
elementary factor sequences over exact field-backed rings, with source refs and
expected rule metadata. The parent optimizer issue #188 remains out of scope.

## Approach Options

Recommended: add one dedicated fixture module and one internal validator, then
register the validator in the internal test group. This follows the existing
catalog pattern, keeps the schema private to tests, and gives later optimizer
issues a stable source of positive and negative examples.

Alternative: extend an existing Park-Woodburn polynomial fixture catalog. That
would mix full-matrix route evidence with factor-sequence rewrite metadata and
make the validator less focused.

Alternative: add production Steinberg helper functions now and test them
against fixtures. That reaches into issue #188 and is explicitly out of scope.

## Chosen Design

Create `test/fixtures/steinberg_optimization_cases.jl` as module
`SteinbergOptimizationFixtureCatalog`. It returns a catalog with `cases` and
`negative_controls`. Each positive case records:

- `id`, `rule_name`, `description`, `source_refs`, and
  `consumer_issue_ids = ("#188",)`;
- `ring_constructor` and `ring` metadata for an ordinary polynomial ring over
  an exact field;
- `factor_metadata` entries containing `(row, col, coefficient)` for each
  elementary factor;
- `factors`, `expected_rewrite_factors`, `original_product`, and
  `rewritten_product`;
- `rewrite_span`, so optimizers can locate the intended local rewrite.

Positive entries cover identity removal, same-position merge, inverse
cancellation as a zero-sum same-position merge, the two nontrivial commutator
relations, and the disjoint commutator identity.

Create `test/internal/steinberg_optimization_fixtures.jl` as the validator. It
loads the fixture module, verifies unique ids, required positive ids, required
rule names, Section 6 source refs, consumer issue metadata, ordinary
field-backed polynomial rings, square matrix factors over one ring, metadata
matching each factor, exact product equality, and negative controls.

Register the validator in `test/runtests.jl` under the internal group.

## Validation Boundary

The validator should accept only test fixtures. It should not export a public
API or call `elementary_factorization`.

It must reject:

- a negative control with mismatched factor rings;
- a negative control whose stored expected product is stale;
- a negative control claiming a commutator rewrite with invalid indices.

The catalog can contain both local positive examples and tempting invalid
examples, but a case is positive only if its original and rewritten factor
products are exactly equal over a single exact field-backed ordinary polynomial
ring.

## Tests

Focused verification commands:

```bash
julia --project=. -e 'include("test/internal/steinberg_optimization_fixtures.jl")'
julia --project=. test/runtests.jl internal
julia --project=. -e 'using Pkg; Pkg.test()'
```

The issue-required command must exit 0. `Pkg.test()` must include the new
internal validator and keep the default public/internal suite green.

## Self-Review

This design is scoped to fixture metadata, validator rules, test registration,
and committed Superpowers documentation. It does not implement a Steinberg
optimizer, export a public API, change `elementary_factorization`, or broaden
Laurent/ToricBuilder support.
