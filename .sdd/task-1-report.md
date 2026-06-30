# Task 1 Report: Evidence-Backed Automatic Quillen Route

## Status

DONE

## Scope Completed

- Replaced the automatic fixture-only polynomial Quillen path with a narrow supplied-evidence route for ordinary-polynomial `3 x 3` elementary matrices.
- Added automatic evidence extraction via `_polynomial_quillen_supplied_evidence_data(A)`.
- Added automatic patch assembly via `_polynomial_quillen_supplied_evidence_patch(A)`.
- Switched automatic route selection and staged-evidence checks to the supplied-evidence route.
- Added expert/public coverage for:
  - fixture-shaped Quillen inputs now returning `QuillenSuppliedEvidencePatchAssembly`;
  - non-fixture multivariate elementary Quillen inputs over `QQ[X,r,g]`;
  - tampered local-sequence and substitution-chain rejection;
  - public unsupported multivariate determinant-one non-elementary rejection.

## RED Evidence

Command:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Observed before implementation:

- Existing automatic Quillen routing still returned `QuillenGlobalPatchAssembly`.
- The new expert assertions failed on `quillen_cert.evidence.quillen_patch isa Suslin.QuillenSuppliedEvidencePatchAssembly`.
- The old patch shape also lacked `base_term_policy`, which confirmed the branch was still on the deterministic fixture path.

## Implementation Notes

- `_polynomial_quillen_elementary_entry(A)` now identifies supported automatic inputs by exact elementary shape.
- `A(0)` is computed by evaluating the selected variable at zero and preserving the remaining generators.
- The automatic local evidence uses the two-open cover over the second generator and its complement.
- Base-term handling is:
  - `:trivial` when `A(0)` is the identity contribution;
  - `:supplied` with a single elementary base-term factor otherwise.
- Route metadata and staged-failure checks were updated to work with `QuillenSuppliedEvidencePatchAssembly`.

## GREEN Evidence

Commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
git diff --check
```

Observed after implementation:

- `test/expert/park_woodburn_route_certificate.jl`: 86/86 passed
- `test/public/factorization_driver_shell.jl`: 61/61 passed
- `test/public/park_woodburn_polynomial_factorization.jl`: 38/38 passed
- `git diff --check`: clean

## Files Changed

- `src/algorithm/factorization.jl`
- `test/expert/park_woodburn_route_certificate.jl`
- `test/public/factorization_driver_shell.jl`
- `test/public/park_woodburn_polynomial_factorization.jl`

## Commit Created

- `382ffe2 feat: gate quillen route on supplied evidence`

## Concerns

- None for Task 1. The implemented boundary remains intentionally narrow: ordinary-polynomial, `3 x 3`, elementary automatic Quillen routing only.

## Review Fix Addendum

### Review fix summary

- Updated `README.md` and `docs/src/index.md` to state the supported #183 boundary precisely: ordinary-polynomial Quillen patching with verified supplied or Murthy-adapter local evidence, exact cover replay, sequence replay, substitution-chain replay, and trivial or supplied `A(0)` base-term handling.
- Removed the stale language that treated #183 as staged and removed any implication of #184, #185, #186, #187, Laurent/ToricBuilder mainline support, or factor-count optimization being available.

### Tests run and outcomes

- `git diff --check` — passed with no whitespace or patch-format issues.

### Commit created

- `docs: clarify quillen #183 boundary`
