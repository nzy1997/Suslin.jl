# Issue 10 Constructive Quillen Patching Design

## Context

Issue #10 depends on the Laurent normalization, Laurent linear solver,
Laurent elementary core, and test-command contract work from issues #8, #9,
#12, and #21. GitHub shows those dependencies are closed as completed, and
the current checkout includes the corresponding APIs:
`normalize_laurent_object`, `lift_laurent_normalization`,
`solve_laurent_linear`, Laurent-aware elementary routines, and the documented
full-suite command `julia --project=. test/runtests.jl all`.

The existing Quillen layer is only scaffolding: `LocalCertificate` stores raw
indices and denominators, `common_denominator_factor` multiplies toy
denominators, and `patched_substitution` performs exact substitution. The new
work should turn that into an inspectable patch object without taking over the
larger factorization driver.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.

## Approaches Considered

1. Add an explicit denominator-coverage patch model. Local data records a
   certificate, denominator, coverage multiplier, and elementary correction.
   Patch construction requires the exact coverage identity
   `sum(multiplier_i * denominator_i) == 1`, builds denominator-weighted
   elementary factors, multiplies them, and verifies the product against the
   requested target correction.
2. Add ad hoc helpers around the current denominator product and substitution
   scaffolding. This would be quick, but it would not produce the inspectable
   data model requested by the issue and would leave verification spread across
   callers.
3. Implement a full Quillen-patching solver that computes coverage
   coefficients internally. This is broader than the issue: denominator
   metadata is already part of the input interface, and internal Bezout
   solving would create new polynomial-only substeps before the patch data
   model is stable.

The chosen design is approach 1. It is the narrowest constructive layer that
uses supplied denominator metadata, exposes stable records for later roadmap
steps, and keeps exact verification local to the patch object.

## Data Model

Extend the Quillen scaffolding with these public records:

- `QuillenDenominatorData(denominator, coverage_multiplier)` records one
  denominator and its supplied coefficient in the coverage identity.
- `QuillenElementaryCorrection(row, col, entry)` records one elementary local
  correction entry. The patcher supports this exact elementary correction
  shape first because Suslin already has verified elementary matrix and factor
  sequence primitives.
- `QuillenLocalContribution(certificate, denominator, coverage_multiplier,
  correction)` combines the local certificate data with the denominator
  metadata and elementary correction.
- `QuillenPatchVerification` records exact verification data:
  denominator-data consistency, coverage sum, product matrix, target matrix,
  and product equality.
- `QuillenPatch` records the ring, matrix size, substitution variable,
  denominator data, local contributions, assembled factors, product, target,
  and verification data.

The existing `LocalCertificate` remains available and is promoted by allowing
abstract integer vectors while continuing to reject mismatched indices and
denominator lengths.

## Construction Semantics

Add:

```julia
construct_quillen_patch(n, X, contributions; target)
```

The constructor:

1. Validates that `target` is an `n x n` matrix over a supported exact
   ordinary polynomial or Laurent polynomial parent.
2. Validates that `X` is one of the target ring generators and records it as
   the substitution variable.
3. Coerces every denominator, coverage multiplier, and elementary correction
   entry into the target ring.
4. Requires every contribution denominator to be present in its
   `LocalCertificate.denominators`.
5. Requires the two elementary indices to appear in
   `LocalCertificate.indices`, because the certificate should witness the
   coordinates used by the local correction.
6. Computes the coverage sum
   `sum(coverage_multiplier_i * denominator_i)`.
7. Throws `ArgumentError("denominator coverage must sum to one")` unless the
   coverage sum is exactly `one(R)`.
8. Builds each exact elementary factor
   `E(row_i, col_i, coverage_multiplier_i * denominator_i * entry_i)`.
9. Multiplies the assembled factor sequence from the identity matrix.
10. Throws `ArgumentError("constructed Quillen patch does not multiply to the
    target correction")` unless the product equals `target`.
11. Returns an inspectable `QuillenPatch`.

For two contributions with the same elementary position and the same local
entry `a`, the product is exact:

```julia
E_ij(c1 * d1 * a) * E_ij(c2 * d2 * a)
    == E_ij((c1 * d1 + c2 * d2) * a)
    == E_ij(a)
```

when `c1 * d1 + c2 * d2 == 1`. This gives the focused constructive patching
example required by the issue without pretending to implement the whole
Quillen-Suslin algorithm.

## Verification Semantics

Add:

```julia
verify_quillen_patch(patch)::Bool
```

Verification recomputes denominator-data consistency, coverage, factor
product, stored product consistency, and product-target equality from the
patch's recorded fields. It returns `false` for invalid or tampered patch
records instead of throwing, except for interrupts.

This supports the negative control where an otherwise valid patch has one
recorded denominator tampered after construction: the stored factor product can
still equal the target, but denominator-data consistency and recomputed
coverage fail.

## Laurent Boundary

Patch assembly itself uses only exact ring coercion and exact matrix
multiplication, so no internal polynomial-only solve is needed. Laurent inputs
therefore remain valid exact-ring inputs. Later polynomial-only Quillen
substeps can normalize Laurent data using the existing issue #8 helpers before
calling this layer and lift back afterward; this patching layer does not hide
that boundary.

## Files

- Modify `src/core/groebner_tools.jl`: promote the certificate scaffolding and
  add denominator data, local contribution, patch construction, and patch
  verification helpers.
- Modify `src/Suslin.jl`: export the new patching data model and helper
  functions.
- Create `test/expert/quillen_patching_exact.jl`: exact constructive patching
  tests for multi-certificate coverage, exact product equality, stable
  denominator tracking, coverage failure, and tampered denominator failure.
- Modify `test/expert/quillen_induction.jl`: keep existing scaffolding tests
  compatible with the promoted certificate constructor.
- Modify `test/runtests.jl`: include the new exact patching expert test.
- Modify `test/public/api_surface.jl`: cover the exported patching API.

## Verification

Issue-specific focused command:

```bash
julia --project=. -e 'include("test/expert/quillen_patching_exact.jl")'
```

Package entry point required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Documented full suite from #21:

```bash
julia --project=. test/runtests.jl all
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design is approved automatically because it selects the conservative explicit
coverage model, keeps the patching layer inspectable, verifies exact
multiplication, supports ordinary polynomial and Laurent parents, and avoids
irreversible API or solver ownership changes.

## Spec Self-Review

- No incomplete markers remain.
- The scope is limited to the constructive denominator-tracking patch object,
  exact assembly, exact verification, and focused tests.
- Negative controls cover missing denominator coverage and tampered
  denominator metadata.
- The design does not broaden into full factorization, determinant
  normalization, or internal Bezout solving.
