# Issue 218 Quillen Supplied Evidence Patch Assembly Design

## Context

Issue #218 is the assembly step after the internal Quillen evidence layers from
#214 through #217. The repository can now represent checked local factor
sequences, extract denominator-cover candidates from those sequences, solve a
powered denominator cover, and replay the Park-Woodburn substitution chain from
`A(X)` to `A(0)`.

The remaining gap is an internal constructor that consumes supplied local
sequence evidence for one ordinary-polynomial matrix and assembles a replayable
global patch without selecting behavior by fixture id. The constructor must
keep the sequence expansion visible: each generated global elementary factor
must be traced to a replayed local sequence factor and a verified cover term,
not to a collapsed opaque local correction.

## Approaches Considered

Recommended: add a new internal supplied-evidence assembly certificate beside
`QuillenGlobalPatchAssembly` in `src/algorithm/quillen_induction.jl`. The
constructor verifies the #214 sequence certificates, extracts the #215
candidate, solves the #216 cover, builds or validates the #217 substitution
chain, expands every local sequence factor using the corresponding solved cover
term, records explicit base-term handling, multiplies the assembled factors,
and verifies exact equality with `A`.

Alternative: retrofit `QuillenGlobalPatchAssembly` so its
`local_certificates` field accepts sequence certificates. This risks weakening
the existing single-correction checks and would force older tests and route
adapter code through a broader union-shaped contract.

Alternative: convert each sequence into one deterministic local realization
certificate and reuse `assemble_deterministic_quillen_patch`. This is too close
to the behavior the issue rejects: the ordered sequence would become one
opaque correction unless the collapse were separately replayed.

The chosen design is the recommended additive internal certificate. It reuses
the existing verifier helpers where they are exact, but adds a sequence-specific
expansion layer so the assembly record exposes the evidence chain required by
the issue.

## Certificate Shape

Add these internal records:

- `QuillenSequenceContributionExpansionVerification`: replay data for one
  verified local sequence expanded under one solved cover term.
- `QuillenSequenceContributionExpansion`: the generated elementary factors for
  one local sequence. It records the source sequence certificate, the solver
  result, local index, powered denominator, coverage multiplier, cover term,
  replayed factor provenance, generated global elementary factors, and replay
  metadata.
- `QuillenSuppliedEvidencePatchAssemblyVerification`: replay data for the full
  constructor output.
- `QuillenSuppliedEvidencePatchAssembly`: the assembled internal patch
  certificate.

The full patch certificate records:

- the original matrix `A`, ring, size, and selected variable;
- the verified #214 local factor sequence certificates;
- the extracted #215 denominator-cover candidate;
- the #216 solver result and cover certificate;
- the #217 substitution chain;
- explicit base-term policy and any supplied base-term factors;
- the per-local sequence expansions;
- the final global elementary factor list, product, target matrix, replay
  metadata, and verification object.

Do not export these names from `src/Suslin.jl`; expert/internal tests use the
`Suslin.` module qualifier.

## Construction And Replay

Add `assemble_quillen_patch_from_local_evidence(A, X, local_certificates; ...)`.
The constructor accepts a nonempty collection of
`QuillenLocalFactorSequenceCertificate` values for the same original matrix,
ring, size, and selected variable. It rejects unverified locals, mixed context,
unproven denominator covers, non-replaying substitution chains, missing
base-term evidence, and final products that do not multiply exactly to `A`.

Construction steps:

1. Validate `A` as a square exact ordinary-polynomial matrix and `X` as a ring
   generator.
2. Verify all local sequence certificates and extract a
   `QuillenDenominatorCoverCandidate`.
3. Solve the candidate with the existing bounded #216 solver. Optional solver
   keyword arguments are passed through so tests can supply wrong multipliers
   and get a deliberate `ArgumentError`.
4. Build a `QuillenPatchSubstitutionChain`, or validate a supplied chain
   against the same matrix, variable, and solver result.
5. Expand each local sequence explicitly. For local index `i`, use the solved
   cover term `g_i * r_i^l`; for each structured factor in that sequence,
   generate `elementary_matrix(n, row, col, cover_term * numerator, R)`. This
   keeps every factor tied to a replayed local sequence factor and the exact
   cover term used for that local sequence.
6. Handle the `A(0)` boundary with an explicit policy:
   `:supplied` requires supplied base-term factors whose product equals the
   chain base term, `:trivial` requires the chain base term to be the identity,
   and `:already_handled` records the staged boundary while requiring the final
   assembled factors to multiply exactly to `A`. Omitting both supplied factors
   and a policy throws a clear `ArgumentError`.
7. Concatenate base-term factors when supplied, then sequence expansion factors,
   multiply exactly, and require the product to equal `A`.

Replay recomputes the candidate, solver cover data, substitution chain,
base-term product, sequence expansions, global factor list, final product, and
metadata from stored inputs. `verify_quillen_patch` gets a method for the new
assembly type that returns `false` on tampered stored data.

## Tests

Create `test/expert/quillen_supplied_evidence_patch_assembly.jl`. The positive
case builds #214 sequence certificates from the supplied local evidence in the
fixture catalog, calls the new constructor with
`base_term_policy = :already_handled`, and checks:

- every child local sequence certificate verifies;
- extracted raw denominators match the local sequence product denominators;
- the solver verifies `sum(g_i * r_i^l) == 1`;
- the substitution chain verifies and telescopes;
- every generated global factor is the expected cover-term expansion of one
  local sequence factor;
- base-term handling is recorded as `:already_handled`;
- the returned factor product and target are exactly `A`;
- the verifier rejects tampering.

Negative controls cover incomplete evidence, corrupted local sequence data,
unproven supplied cover multipliers, tampered substitution chains, and omitted
base-term policy/evidence.

Update `test/expert/quillen_induction_constructive.jl` with a focused
regression that the constructive acceptance path can also use the supplied
evidence constructor, without replacing the existing deterministic assembly
tests.

Register the new expert test in `test/runtests.jl`.

## Files

- Modify `src/algorithm/quillen_induction.jl`.
- Create `test/expert/quillen_supplied_evidence_patch_assembly.jl`.
- Modify `test/expert/quillen_induction_constructive.jl`.
- Modify `test/runtests.jl`.
- Add `docs/superpowers/plans/2026-06-30-issue-218-quillen-supplied-evidence-patch-assembly.md`.

No Murthy local solving, general `SL_3` driver, ECP implementation, new public
API export, route-boundary wiring, or invented base-term factorization
algorithm is included.

## Verification

Focused supplied-evidence assembly:

```bash
julia --project=. -e 'include("test/expert/quillen_supplied_evidence_patch_assembly.jl")'
```

Existing constructive Quillen checks:

```bash
julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'
```

Package entry point:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Diff hygiene:

```bash
git diff --check
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
recommended additive internal certificate is approved automatically because it
follows the issue's supplied-evidence objective, keeps sequence expansion
explicit, preserves existing public APIs, and leaves route-boundary work to
#220.

## Spec Self-Review

- No placeholders remain.
- The design covers every required replay stage from local sequence
  verification through exact final product equality.
- Base-term evidence is explicit and deliberately staged when not supplied.
- Sequence factors are generated from structured factor fields and solved cover
  terms, not from fixture ids or opaque corrections.
- Negative controls cover all staged errors named in the issue.
