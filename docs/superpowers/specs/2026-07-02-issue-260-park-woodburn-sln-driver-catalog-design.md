# Issue 260 Park-Woodburn SLn Driver Catalog Design

Issue #260 adds a reusable ordinary-polynomial `SL_n` recursive driver fixture
catalog for parent #186. The catalog records source-grounded matrices and
metadata for the intended chain from last-column ECP reduction, through
right-clearing and recursive descent, to a final evidence-backed `SL_3` route.
It does not implement recursive `SL_n` factorization or broaden public
Park-Woodburn support.

## Context

There is no repository `AGENTS.md`; the README test instructions apply. The
current branch starts from `main` after the #184 `SL_3` driver and #185 ECP
acceptance work. Existing catalogs follow a two-file pattern: a fixture module
under `test/fixtures/` and an internal validator under `test/internal/`, with
`test/runtests.jl` including the validator in the `internal` group.

GitHub issue/PR fetching is unavailable in this Agent Desk sandbox because the
GitHub CLI cannot reach the configured proxy. The issue body supplied by the
run, local merge history, README scope text, and existing fixture validators are
the source of truth for this design.

## Chosen Approach

Create a new catalog instead of extending
`test/fixtures/park_woodburn_polynomial_cases.jl`.

1. Add `test/fixtures/park_woodburn_sln_driver_cases.jl` with explicit #186
   support roles: `:issue186_mainline`, `:staged_issue186_candidate`, and
   `:legacy_regression`.
2. Add `test/internal/park_woodburn_sln_driver_fixtures.jl` to validate the
   catalog and negative controls without relying on diagnostic string parsing.
3. Register the validator in `test/runtests.jl`.

This is preferable to expanding the legacy recursive column-peel catalog
because #186 needs route-provenance, peel-step, ECP, right-clearing, final
`SL_3`, staged-boundary, and parent-gate metadata that would make the older
fast-local/disjoint-block route schema ambiguous.

## Catalog Surface

Each positive or staged-positive entry records:

- A unique id and ordinary exact field-backed polynomial matrix of size
  `n >= 4`.
- `support_role` and `route_provenance`, so a legacy fast-local or
  disjoint-block regression cannot be mistaken for #186 mainline support.
- `expected_peel_count`, `descent_dimensions`, and per-step metadata.
- Last-column ECP metadata, including whether the column evidence is replayed,
  staged, or legacy-only.
- Right-clearing metadata that replays the expected `after_left_matrix`,
  `right_factors`, `peeled_matrix`, and `next_block`.
- Final `SL_3` route metadata, distinguishing evidence-backed #184 routes from
  legacy shortcuts and missing evidence.
- `staged_reason_codes` as tuples of symbols from a closed set:
  `:missing_ecp_evidence`, `:missing_final_sl3_route`,
  `:unsupported_coefficient_ring`, and `:legacy_regression_only`.
- `source_refs` including
  `refs/arXiv-alg-geom9405003v1 Section "Reduction to SL_3(k[x_1,...,x_m])"`,
  `refs/arXiv-alg-geom9405003v1 Section 4`, and
  `refs/arXiv-alg-geom9405003v1 Section 5` where applicable.
- `consumer_issue_ids = ("#186",)`.

The catalog will include at least five positive or staged-positive cases:

- A legacy recursive column-peel regression preserved from the older catalog
  and marked `:legacy_regression`.
- A multivariate `SL_4` mainline case whose last column uses #185 ECP evidence.
- A multistep `SL_5` case that peels twice before `SL_3`.
- A case whose final block is accepted by the #184 evidence-backed `SL_3`
  route.
- A staged candidate whose ECP or final `SL_3` evidence is explicitly missing.

## Validator

The validator checks the catalog as data, not as an implementation of #186.
It rejects:

- Duplicate ids across cases and negative controls.
- Non-polynomial rings, inexact rings, or non-field-backed coefficient rings.
- Positive/staged-positive matrices with size `< 4`.
- Positive/staged-positive matrices whose determinant is not one.
- Inconsistent peel counts, descent dimensions, or step dimensions.
- Last-column metadata whose column does not match the step input matrix.
- Replayed ECP metadata without a verified `ECPColumnReductionCertificate`.
- Right-clearing metadata that does not replay the stored factors and next
  block exactly.
- Final `SL_3` support claims without evidence-backed #184 route metadata.
- Unknown or malformed staged reason codes.
- Entries claiming `:issue186_mainline` support while using only legacy
  fast-local/disjoint-block provenance.

Negative controls cover determinant-not-one input, unsupported coefficient
ring, corrupted peel expectation, unknown staged reason code, and a false #186
support claim without required ECP or final `SL_3` evidence.

## Tests

The target verification command is:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_sln_driver_fixtures.jl")'
```

The package gate remains:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The internal validator should also be included by the default test runner
through `test/runtests.jl`.

## Non-Claims

This issue does not implement ECP, `SL_3` realization, recursive peeling,
public dispatch, Laurent/ToricBuilder support, #187 final acceptance, or
Steinberg factor-count optimization. Staged entries are reusable planning and
regression inputs, not acceptance claims.

## Self-Review

The design is focused on a new fixture catalog and validator. It preserves the
legacy recursive column-peel cases as regression references while requiring
explicit role/provenance metadata before any case can count as #186 mainline
support. The negative controls directly match the issue's requested rejection
cases.
