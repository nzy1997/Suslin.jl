# Issue 270 Park-Woodburn Mainline Acceptance Catalog Design

Issue #270 adds the final public acceptance catalog for parent #187. The
catalog records the ordinary-polynomial Park-Woodburn mainline boundary without
changing public factorization behavior or treating older regression fixtures as
the full public claim.

## Context

There is no repository `AGENTS.md`; the README and existing test catalog
patterns apply. Live GitHub issue fetching is unavailable in this Agent Desk
sandbox because `gh issue view` cannot reach the configured proxy, so the
supplied #270 issue body, local merge history through #266, and the existing
fixture validators are the source of truth.

The repository now has separate evidence catalogs for the major #187 layers:

- #181/#182 normality and Murthy local `SL_3` support boundaries.
- #183 Quillen patching mainline fixtures.
- #184 public `SL_3` driver cases in
  `test/fixtures/park_woodburn_sl3_driver_cases.jl`.
- #185 ECP mainline cases in `test/fixtures/ecp_mainline_cases.jl`.
- #186 recursive `SL_n`, `n > 3`, driver cases in
  `test/fixtures/park_woodburn_sln_driver_cases.jl`.

## Chosen Approach

Use a new final catalog and validator:

1. Create `test/fixtures/park_woodburn_mainline_acceptance_cases.jl` with
   concise end-user `elementary_factorization(A)` cases. Entries record their
   public route expectation, expected status, source refs, upstream issue
   coverage, and linked upstream evidence ids.
2. Create `test/internal/park_woodburn_mainline_acceptance_fixtures.jl` as the
   validator and smoke test. It validates the new schema, positive case ring
   and determinant conditions, source references to Park-Woodburn Sections 2-5,
   upstream issue coverage, route/status metadata, and negative controls.
3. Register the internal validator in `test/runtests.jl` so `Pkg.test()` runs
   it with the rest of the internal fixture gates.

This is preferable to extending `park_woodburn_polynomial_cases.jl` because
that older catalog remains useful regression input but does not separate legacy
fixture support from the full #187 mainline claim. It is also preferable to
changing production code because #270 is an acceptance catalog and explicitly
leaves ECP implementation, recursive implementation, Laurent/ToricBuilder
support, and factor-count optimization out of scope.

## Catalog Shape

Each case uses a NamedTuple schema with:

- `id`, `entry_class`, `expected_status`, and `public_route`.
- `ring_constructor`, `ring`, `matrix`, and `determinant_metadata`.
- `source_refs`, including the Park-Woodburn section exercised by the case.
- `upstream_issue_ids` and `upstream_evidence` metadata.
- `acceptance_metadata` explaining whether the case is a supported public
  example, staged/missing-evidence boundary, or negative control.

The required entry classes are:

- an evidence-backed multivariate `SL_3` case consuming #184;
- an `SL_n`, `n > 3`, case consuming #185 ECP and #186 recursive evidence;
- a README-style ordinary-polynomial example;
- a determinant-not-one negative control;
- an unsupported coefficient-ring negative control;
- a missing-evidence negative control.

Positive supported cases must use exact field-backed ordinary polynomial rings,
have determinant one, expose source refs for the relevant Park-Woodburn layers,
and include the upstream evidence ids required by their selected route. Staged
cases may be determinant-one but must name the missing evidence explicitly and
must not claim #187 support.

## Validation Boundary

The validator should reject:

- duplicate ids across cases and negative controls;
- positive entries over non-field-backed or non-exact ordinary polynomial rings;
- positive entries whose determinant metadata is false or whose determinant is
  not one;
- `:mainline_accepted` route expectations that omit required #184, #185, or
  #186 upstream evidence ids;
- entries whose source refs omit the Park-Woodburn section for the layer being
  exercised;
- negative controls that unexpectedly validate.

The validator should not call `elementary_factorization`; it is a catalog
acceptance gate. Existing public tests continue to cover behavior.

## Tests

Focused verification commands:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_mainline_acceptance_fixtures.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

The first command is the issue-required validator check. The full package test
ensures the new internal fixture is registered and does not disturb existing
public/internal acceptance tests.

## Self-Review

The design is scoped to fixture metadata, validation, and test registration.
It does not add public APIs, close #187 by itself, change factorization
behavior, broaden coefficient-ring support, claim Laurent/ToricBuilder support,
or optimize factor counts. The staged boundary is explicit when recursive or
other upstream evidence is missing.
