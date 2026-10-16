# Issue 355 Park LaurentToPoly Conversion Design

## Scope

Issue #355 now ends at a verified Laurent-to-polynomial conversion certificate.
It must not construct monicity plans, transformed polynomial columns, ordinary
ECP children, lifted ordinary factors, production dispatch, or full d14 replay
to `e14`.

## Chosen Approach

Use the generic Park Algorithm 6.1 preprocessing path already prototyped in the
worker branch:

1. Preserve existing already-polynomial and Laurent-unit fast paths.
2. For non-unit, non-polynomial two-generator Noether columns, normalize the
   lowest selected-variable coefficient to one by determinant-one elementary
   factors.
3. Clear nonpositive selected-variable degrees from every other entry by
   elementary row additions.
4. Recompute the monomial polynomialization power from the transformed support.
5. Build the ordinary polynomial column and inverse lift maps.
6. Store a polynomial Bezout witness and verify `sum(w[i] * column[i]) == 1`
   directly during certificate replay.

The d14 fixture selection `selected_entry = 9`, `selected_generator = v` remains
in the report/test layer. The generic code receives only the selected entry and
generator supplied by the caller.

## Alternatives Considered

The old full-bridge route would continue into monicity search and ordinary ECP
construction. That is now #366-#369 and is out of scope.

A production reducer fallback is also out of scope. This PR keeps the existing
public dispatch behavior unchanged and only extends the internal conversion
certificate.

## Verification Design

Internal tests cover all three conversion families: already polynomial, Laurent
unit entry, and generic Park preprocessing without a unit. Negative controls
tamper with preprocessing factors, intermediate columns, polynomialization
power, maps, polynomial columns, and witness coefficients.

The bounded d14 report constructs the checked-in `case_008` d14 column, runs the
generic conversion with entry 9 and generator `v`, verifies the certificate, and
records finite source/polynomial/witness term metrics and elapsed time.
