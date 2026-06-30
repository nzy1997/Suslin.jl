# Issue 220 Quillen Patching Gate Design

## Context

Issue #220 closes the #183 Quillen local-to-global patching parent. The current
checkout already contains the dependency work through #219: sequence
certificates, denominator extraction, cover solving, Park-Woodburn substitution
chains, supplied-evidence patch assembly, and Murthy-adapter consumption.

The remaining gap is the ordinary-polynomial factorization boundary. It still
has an automatic Quillen shortcut that recognizes one hard-coded fixture-shaped
matrix and builds a deterministic patch directly from local realization
certificates. That path proves less than the closeout gate requires because it
does not force route acceptance to pass through replayable local sequence
evidence, exact cover solving, substitution-chain replay, and base-term policy
handling.

## Selected Approach

Use the existing supplied-evidence assembler as the automatic Quillen route for
the supported ordinary-polynomial boundary. The automatic route will build
verified `QuillenLocalFactorSequenceCertificate` inputs, call
`assemble_quillen_patch_from_local_evidence`, and then adapt the resulting
`QuillenSuppliedEvidencePatchAssembly` through the existing
`PolynomialQuillenPatchRouteAdapter`.

The supported automatic boundary is intentionally narrow:

- exact ordinary polynomial rings over supported coefficient fields;
- `3 x 3` elementary matrices over an ordinary polynomial ring with at least two
  generators;
- the first generator is the Park-Woodburn substitution variable;
- the second generator supplies the two-open cover `s, 1-s`;
- local evidence covers the non-base delta `A / A(0)`;
- `A(0)` is accepted only when it is trivial or supplied as explicit elementary
  base-term factors by the route builder.

This retains #183 acceptance without implementing #184 general `SL_3`, #185
ECP, #186 recursive `SL_n`, #187 full public acceptance, Laurent/ToricBuilder
mainline acceptance, or factor-count optimization.

## Data Flow

Automatic polynomial factorization tries existing fast-local and block routes
first. If those do not apply, it asks the Quillen supplied-evidence builder for
a patch.

For a supported elementary matrix `A = E_ij(f)`, the builder computes
`base = A(0)` by substituting the selected variable with zero and `delta` by
subtracting the base entry from `f`. It creates two local sequence certificates
for denominators `s` and `1 - s`. Each certificate has one replayable elementary
local factor for the delta term and records provenance showing it came from the
automatic #183 supplied-evidence boundary, not from a fixture id.

The builder handles `A(0)` through the supplied-evidence base-term policy:

- if `base` is the identity, it uses `base_term_policy = :trivial`;
- if `base` is a single elementary matrix, it supplies that matrix as
  `base_term_factors` with `base_term_policy = :supplied`;
- otherwise it rejects the route with a staged message mentioning missing
  `A(0)` base-term evidence.

The produced patch must verify with `verify_quillen_patch(patch)` and its global
elementary factors must multiply exactly to the route target before
`elementary_factorization` returns them.

## Error Boundary

Unsupported determinant-one multivariate matrices continue to fail before
returning factors. The failure message must state that verified Quillen/local
evidence is missing. Unsupported coefficient rings and determinant-not-one
inputs keep the existing earlier failures.

Tampered local sequence certificates, denominator covers, substitution chains,
base-term factors, supplied/Murthy adapter records, and Quillen patch
certificates must fail verification or assembly before route certificates can
return factors.

## Tests And Documentation

Focused expert tests will add #220 acceptance assertions that the automatic
Quillen route uses a `QuillenSuppliedEvidencePatchAssembly`, sequence
certificates replay factor-by-factor, cover solving replays from raw
denominators through powered terms, substitution chains telescope and record
`A(0)`, and tampering fails.

Public tests will include a non-fixture elementary Quillen case accepted through
the new supplied-evidence route and verify
`verify_factorization(A, factors) == true`. Negative controls will cover
determinant-not-one input, missing local evidence, unsupported coefficient ring,
tampered sequence data, tampered substitution chains, missing required `A(0)`
evidence, and tampered patch certificates.

README and docs wording will state the precise supported #183 boundary:
ordinary-polynomial Quillen patching with verified supplied or Murthy-adapter
local evidence, exact cover replay, sequence replay, substitution-chain replay,
and trivial or supplied base-term handling. The docs will explicitly keep
#184-#187, Laurent/ToricBuilder mainline support, and factor-count optimization
staged.
