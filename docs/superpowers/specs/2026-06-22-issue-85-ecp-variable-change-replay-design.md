# Issue 85 ECP Variable-Change Replay Design

## Context

Issue 84 added a replayable internal certificate path for ECP column
reductions. Its monicity-normalization route already performs the bounded
Park-Woodburn-style substitution search, reduces the substituted column through
the existing supported unit/witness stages, inverse-substitutes those factors,
and verifies that the lifted factors reduce the original column to `e_n`.

Issue 85 narrows the next step: record the variable-change data explicitly
enough for exact replay. The issue comment is binding: this must not become a
new search algorithm, and the record must state whether the selected monic entry
has been moved or normalized into the first coordinate required by the later
Park-Woodburn link theorem.

There are no repository instruction files such as `AGENTS.md`, `CLAUDE.md`, or
`GEMINI.md` in this worktree. The README says the package does not commit a
`Manifest.toml`, and the baseline `Pkg.test()` passed on the unmodified branch.

## Design Choice

Extend the existing `:monicity_normalization` stage in
`src/algorithm/column_reduction.jl` so it is the explicit variable-change stage
record. Keep the current public behavior:

- `reduce_unimodular_column(v, R)` still returns only factors.
- `ecp_column_reduction_certificate(v, R)` remains the expert certificate path.
- No public export is added.
- No general Quillen patching or broader variable-change search is added.

Alternatives considered:

- Add a separate variable-change reducer. This would duplicate the current
  deterministic bounded substitution path and risk factor/certificate drift.
- Change the public reducer to return a rich stage object. This violates the
  issue's public-entry-point constraint.
- Rename the stage kind away from `:monicity_normalization`. Keeping the
  existing kind avoids churn in the issue 84 certificate tests while the added
  fields make the variable-change data explicit.

## Stage Record

For the deterministic supported stage, the record stores and replay checks:

- the original input column;
- the variable order as the ring generators;
- the selected source variable/index, target variable/index, shift sign, shift
  power, and shift polynomial;
- forward and inverse substitution maps, plus their ordered value vectors;
- the transformed column;
- the selected monic entry index and value in the transformed column;
- the transformed substage and transformed-stage factors;
- inverse-substituted factors over the original ring;
- final factors and output column over the original ring;
- a first-coordinate status that is `:already_first` when the selected monic
  entry is the first coordinate and `:not_moved` otherwise, with an explicit
  empty move-factor list for the current implementation;
- verification fields for the selected monic check, transformed reduction, and
  original reduction.

Replay recomputes all of these values from the recorded scalar choices. Extra,
missing, or corrupted fields must make `verify_ecp_column_reduction(cert)`
return `false`.

## Fixture Scope

The issue requires at least two catalog cases that avoid direct unit entries
and the current witness-unit shortcut until variable change creates a monic
entry. The merged catalog exposes one such route. Add one narrow second catalog
case by permuting the existing GF(2) variable-change example, so the algebra is
already covered but the stage replay sees a distinct catalog id and selected
monic coordinate.

The new catalog entry remains test fixture data. The reducer behavior is not
broadened.

## Tests

Add `test/expert/ecp_variable_change_replay.jl` and register it in the expert
test group. The focused test will:

- build certificates for `ecp-variable-change-monic-gf2` and the new permuted
  variable-change catalog case;
- assert the original columns have no unit entry and the current direct
  supported unit/witness shortcut returns `nothing`;
- assert the recorded forward substitution rebuilds the transformed column;
- assert the recorded selected monic index identifies a monic transformed
  entry;
- assert transformed-stage factors reduce the transformed column to `e_n`;
- assert inverse-substituted factors are over the original ring and reduce the
  original column to `e_n`;
- assert the stage verifier accepts the unmodified certificate;
- corrupt the inverse substitution map and selected monic index and assert
  replay rejects the certificate.

Existing `test/expert/ecp_column_certificate.jl` remains the broad certificate
coverage. Its monicity tamper helper is updated if field names or expected keys
change.

## Verification

Focused issue command:

```bash
julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expert regression command:

```bash
julia --project=. test/runtests.jl expert
```

## Spec Self-Review

- No incomplete markers remain.
- Scope stays inside the existing certificate/replay route.
- The public reducer return type is unchanged.
- Every newly recorded field has a replay check.
- The first-coordinate Park-Woodburn guardrail is recorded but no link-theorem
  move is implemented in this issue.
