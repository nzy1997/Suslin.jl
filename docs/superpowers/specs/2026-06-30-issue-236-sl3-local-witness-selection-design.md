# Issue 236 SL3 Local Witness Selection Design

## Context

Issue #236 builds on the checked #235 `SL3RealizationInputContext`. The new
layer must choose a selected induction variable and the local Murthy
special-form evidence for ordinary-polynomial `SL_3` driver contexts. It must
not claim that arbitrary multivariate determinant-one matrices can be moved to
Park-Woodburn local form without replayable evidence.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
This worktree is already isolated on
`agent/issue-236-select-monic-variable-and-local-form-witnesses-f-run-1`.
`gh issue view` for #236/#235/#234 is blocked by the sandbox proxy, so the
Agent Desk issue body plus the landed #234/#235 specs and PR #250 history are
the authoritative context for this run. After `Pkg.instantiate()`, the #235
focused context test passes with 60/60 assertions.

## Approaches Considered

1. Add an internal witness-selection record with constructor and verifier in
   `src/algorithm/factorization.jl`, plus a small reusable monicity-proof
   helper in `src/algorithm/sl3_local.jl`. This is the chosen approach because
   it is additive, internal, replay-focused, and keeps the public route
   dispatcher unchanged.
2. Fold witness selection into `SL3RealizationInputContext`. This would reduce
   the number of records, but it would blur #235's input validation boundary
   with #236's local Murthy witness boundary.
3. Start an automatic coordinate-change search for arbitrary multivariate
   `SL_3` matrices. This is out of scope because #236 requires replayable
   evidence and explicitly forbids broad unsupported searches.

## Design

Add non-exported production helpers:

```julia
Suslin._select_sl3_local_form_witness(context; selected_variable = nothing,
    local_form_witness = nothing, variable_change_metadata = nothing,
    normality_conjugation_metadata = nothing)

Suslin._verify_sl3_local_form_witness_selection(record)::Bool
```

The record stores the #235 context, selected variable, selected-variable index
and name, the extracted special-form entries `(p, q, r, s)`, the corresponding
local-form matrix, a monicity witness for `p` in the selected variable, the
local-form witness metadata used, supplied variable-change metadata, supplied
normality/conjugation metadata, per-source replay statuses, an overall
`support_status`, a `witness_source`, and a staged diagnostic.

Selection is conservative:

- the #235 context must verify before the selector accepts it;
- a selected-variable hint may fill a missing context variable, but any hint
  must be a generator of the context ring and must agree with an already stored
  context variable and index;
- already-special-form matrices are accepted by extracting
  `_sl3_local_target_entries(context.matrix)` and proving that the extracted
  `p` is monic in the selected variable;
- explicit local-form witness metadata may restate `p`, `q`, `r`, `s`,
  `entries`, or a `local_form_matrix`; if supplied, every recorded entry must
  match the extracted local form;
- variable-change or normality/conjugation metadata is replayable only when it
  has a replay payload, a source/input matrix equal to the context matrix, an
  aligned selected variable, and a supplied local-form matrix or local-form
  entries whose `p` is monic in the selected variable;
- metadata shells without a replay payload or without local-form replay data
  are recorded as staged diagnostics, not support.

The monicity proof is a replayable named tuple computed from the current ring:
`(variable, variable_index, degree, leading_coefficient, is_monic)`. It is
verified by recomputing the leading coefficient with the existing
`_sl3_local_coefficient_in_variable_degree` helper.

## Error Handling

The selector throws `ArgumentError` for malformed or contradictory supplied
evidence: non-generator selected-variable hints, selected-variable conflicts,
local-form entries that do not reconstruct the recorded local-form matrix,
local-form `p` that is not monic in the selected variable, source matrices that
do not match the #235 context, variable metadata that names a different
selected generator, and replay metadata that claims a local-form transform but
does not supply enough data to verify it.

When evidence is simply absent, the selector returns a staged record. The staged
diagnostic uses the stable reason `missing supported local-form witness` and
records missing/partial sources so later driver code can stop before consuming
Murthy route inputs.

## Verification

Create `test/expert/park_woodburn_sl3_witness_selection.jl` and register it in
the expert group in `test/runtests.jl`. Focused tests cover:

- the #234 univariate fast-local case as an already-special-form witness;
- a multivariate special-form case with `p` monic in the selected variable;
- exact replay of `(p, q, r, s)`, selected variable metadata, and monicity
  witness fields;
- a supplied variable-change replay witness that aligns with the #235 context
  and supplies a replayable local-form target;
- a supplied normality/conjugation shell without local-form replay data, which
  remains staged;
- a determinant-one `SL_3` matrix with no supported local-form witness, which
  remains staged with the stable diagnostic;
- negative controls for a non-generator selected variable, non-monic `p`,
  corrupted local-form entry metadata, corrupted variable-change metadata, and
  removed required replay payload.

Required commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_witness_selection.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
recommended conservative design is approved automatically because it is
internal, verifier-led, accepts only already-special-form or supplied replayable
local-form data, and stages all unsupported determinant-one contexts.

## Spec Self-Review

- No placeholder markers remain.
- The design separates #235 input validation from #236 witness selection.
- The supported/staged boundary is explicit and does not trust fixture ids or
  arbitrary multivariate matrices.
- Every issue verification case and negative control is represented.
- Public dispatch, Murthy certificate construction, Quillen patching, broad
  coordinate-change search, ECP, and recursive `SL_n` remain out of scope.
