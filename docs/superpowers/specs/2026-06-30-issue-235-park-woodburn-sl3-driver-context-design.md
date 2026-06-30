# Issue 235 Park-Woodburn SL3 Driver Context Design

## Context

Issue #235 adds the checked internal input context that parent #184 needs
before the ordinary-polynomial `SL_3` driver starts choosing variables or
requesting Park-Woodburn evidence. Issue #234 already added a source-grounded
catalog with univariate fast-local, multivariate local-form, Quillen-mainline,
legacy patched-substitution, and staged determinant-one cases. The current
factorization route still checks size, ordinary polynomial rings, determinant
one, selected variables, and staged boundaries in scattered helper functions.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
The GitHub issue fetch through `gh issue view` is blocked by the sandbox proxy,
so the Agent Desk issue body is the authoritative issue text for this run.
Baseline `julia --project=. -e 'using Pkg; Pkg.test()'` passed before changes
with public 382/382 and internal 641/641.

## Approaches Considered

1. Add an internal context record in `src/algorithm/factorization.jl`, with a
   constructor and verifier that recompute all stored fields. This is the
   chosen approach because it centralizes the #184 boundary without changing
   public dispatch or implementing new Park-Woodburn proof layers.
2. Extend the #234 fixture validator only. This would validate catalog data,
   but it would leave production code without a stable checked context for
   later #236-#238 driver layers.
3. Route `elementary_factorization` through the context immediately. This is
   out of scope because the issue explicitly says public behavior must not
   change yet.

## Design

Add non-exported production helpers in `src/algorithm/factorization.jl`:

```julia
Suslin._sl3_realization_input_context(A; selected_variable = nothing,
    catalog_metadata = (;), local_form_witness = nothing,
    variable_change_metadata = nothing,
    normality_conjugation_metadata = nothing,
    quillen_murthy_metadata = nothing)

Suslin._verify_sl3_realization_input_context(context)::Bool
```

The context record stores the input matrix, base ring, coefficient ring, size,
ring profile, generator list and names, selected variable and index, determinant
and determinant status, exact-field status, optional catalog metadata, local
form witness, optional replay metadata for variable-change, normality or
conjugation, and Quillen/Murthy evidence, evidence availability status,
support/staged status, and a replayable staged diagnostic.

Construction deliberately rejects inputs that are not `3 x 3`, ordinary
polynomial, exact field-backed, and determinant-one. It uses the existing
`_validate_factorization_matrix`, `_factorization_ring_profile`, and
`_require_polynomial_sl_determinant` helpers. If a selected variable is
supplied, it must be one of `gens(R)`; #234 selected-variable named tuples are
accepted by reading their `generator` and `index` fields. For univariate
fast-local inputs, the generator can be inferred.

Evidence classification records availability without manufacturing missing
proofs:

- `local_form_status = :fast_local` for the existing univariate fast-local
  shape, or `:replayed` for a supplied local-form witness whose matrix entry is
  monic in the selected variable;
- generic variable-change, normality/conjugation, and Quillen/Murthy metadata
  are `:replayed` only when they include a nonempty replay identifier and a
  target matrix equal to the input;
- metadata with identifiers but no replay target is `:recorded`, which is
  useful diagnostics but does not make the context supported;
- absent evidence is `:missing`.

The context is `:supported` only when at least one evidence path is replayable
or the input is in the existing univariate fast-local slice. Otherwise it is
`:staged`, with a diagnostic naming missing evidence classes and any partial
recorded metadata. If `catalog_metadata.expected_status` is supplied, it must
match the computed support/staged status; a catalog entry cannot mark an
arbitrary determinant-one multivariate `SL_3` input as supported without
replayable evidence metadata.

## Verification

The verifier recomputes every stored field from the matrix and stored hints:
size, base ring, coefficient ring, ring profile, generators, selected-variable
membership in `gens(R)`, determinant one, exact-field status, evidence
availability, support/staged status, and staged diagnostics. It catches ordinary
exceptions and returns `false`, matching existing internal verifier style.

Focused tests in `test/expert/park_woodburn_sl3_driver_context.jl` cover:

- the #234 univariate fast-local case;
- a multivariate `SL_3` case with only recorded Quillen/Murthy ids, which stays
  staged until replay metadata is available;
- a witness-backed multivariate local-form case that is supported;
- a determinant-one multivariate case with no replayable evidence, which stays
  staged;
- constructor rejection for determinant not one, unsupported coefficient ring,
  and non-generator selected variable;
- verifier rejection after corrupting selected variable, determinant status,
  ring profile, evidence status, and staged diagnostics.

Register the expert test in `test/runtests.jl`. Do not export the context type
or route public `elementary_factorization` through it.

## Verification Commands

Focused command required by the issue:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_driver_context.jl")'
```

Package command required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Diff hygiene:

```bash
git diff --check
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it is additive, internal, verifier-led,
keeps public behavior unchanged, and chooses the conservative staged boundary
whenever replayable evidence is absent.

## Spec Self-Review

- No incomplete markers remain.
- The supported/staged boundary is explicit and does not trust fixture ids
  alone.
- Every issue verification case and negative control is represented.
- The scope excludes local witness selection, Murthy calls, Quillen patching,
  coordinate-change search, normality/conjugation proof, and public dispatch.
