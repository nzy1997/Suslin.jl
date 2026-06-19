# Issue 13 Exact Unimodular Column Reduction Design

## Context

`reduce_unimodular_column(v, R)` currently validates that `v` is unimodular and
then tries two small constructive paths: a direct unit entry or a witness row
with a unit coefficient. If those fail, it runs a bounded monicity substitution
search. That keeps earlier examples working, but the public routine does not
yet expose a clear exact pipeline, does not handle longer supported columns
well, and does not use the Laurent normalization and solver decisions from the
dependency issues.

Issue #13 depends on the Laurent normalization helpers from #8, the
Suslin-owned Laurent linear solving decision from #9, Laurent elementary core
support from #12, ToricBuilder fixture constraints from #19, and the test
command contract from #21. The implementation must stay narrower than general
Laurent factorization.

## Approaches Considered

1. Add a staged exact pipeline around the existing constructive kernels. This
   is the chosen approach. It keeps the public function name, makes validation
   and unsupported failures explicit, supports longer columns when a reducible
   3-entry block is present, and handles Laurent columns by normalizing to the
   existing ordinary polynomial helper layer and lifting the exact factors back.
2. Implement a broad Quillen-Suslin or Laurent factorization algorithm now.
   This would exceed the issue scope and overlap later roadmap issues that own
   full matrix factorization and general Laurent acceptance.
3. Keep the heuristic behavior and add only tests. That would preserve old
   behavior but would not deliver the requested exact length 6, 8, and 12
   reductions or the normalized Laurent path.

## Design

`reduce_unimodular_column(v, R)` remains the public entry point. It will be
structured internally as these stages:

- Input validation and coercion: require one-based indexing, length at least
  three, coercibility into `R`, and immediate `ArgumentError` if
  `is_unimodular_column(column, R)` is false.
- Ring profile selection: ordinary polynomial rings use the polynomial
  pipeline; Laurent polynomial rings first try direct Laurent unit reduction and
  then use `normalize_laurent_object(column)`.
- Witness extraction: keep `_unimodular_witness(column, R)` as the exact
  Groebner/ideal witness provider and use it only inside stages that can turn
  the witness into constructive elementary factors.
- Exact reduction steps: preserve the direct unit and witness-unit paths, keep
  the existing finite monicity substitution path as a checked exact small-column
  stage, and add a 3-entry block stage for longer supported columns. The block
  stage scans 3-entry subcolumns, reduces a supported unimodular subcolumn to
  the local `e_3`, embeds those factors into the full dimension, and then kills
  all outside entries using the created pivot unit.
- Laurent normalization: for Laurent columns without an immediate unit
  reduction, multiply by the monomial shift recorded by
  `normalize_laurent_object`, reduce the ordinary polynomial column exactly,
  lift each factor entry back into the Laurent ring, and prepend unit
  normalization factors that convert the resulting monomial-unit last entry to
  one.
- Staged unsupported failure: if validation succeeds but no exact stage applies,
  throw an `ArgumentError` whose message identifies unsupported exact
  unimodular-column reduction and does not claim the input is non-unimodular.

Every stage that returns factors must be checked by exact multiplication before
the public function returns them. A failed internal check is an implementation
error, not a user-facing unsupported case.

## Tests

Add `test/expert/unimodular_reduction_exact.jl` and register it in the expert
test group. The file covers:

- Supported ordinary polynomial columns of lengths 6, 8, and 12 reducing
  exactly to `e_n`.
- A normalized Laurent column with negative exponents reducing exactly to
  `e_n`.
- The existing small examples from `test/expert/unimodular_columns.jl` still
  reducing exactly.
- A non-unimodular column failing immediately with `ArgumentError`.
- A unimodular but unsupported column failing with a staged unsupported message
  that does not say the column is non-unimodular.

The issue-specific command is:

```bash
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
```

The required Agent Desk package entry point is:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The documented full suite from #21 is:

```bash
julia --project=. test/runtests.jl all
```

## Scope Boundaries

This design does not claim general Laurent factorization, does not change the
public function name, does not add new public API, and does not introduce a new
solver backend. The Laurent path reuses existing normalization metadata and
ordinary polynomial reduction stages; unsupported cases remain explicit.

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design is approved automatically because it chooses the narrowest staged exact
pipeline that satisfies issue #13, preserves existing small behavior, reuses
the dependency issue helpers, and avoids irreversible public API changes.

## Spec Self-Review

- No incomplete markers remain.
- The scope is focused on `reduce_unimodular_column` and expert coverage.
- Unsupported unimodular inputs are distinct from non-unimodular validation
  failures.
- Laurent handling is explicitly normalized and lifted back, without claiming
  broader Laurent matrix factorization support.
