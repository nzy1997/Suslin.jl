# Issue 244 ECP Monicity Normalization Design

## Context

Issue #243 added checked ordinary-polynomial ECP input contexts. Issue #244 adds the next replayable stage: pick a selected variable, make a selected column entry monic in that variable by either identity or bounded deterministic variable change, and move the selected monic entry into coordinate 1 when needed.

The current reducer already has `_deterministic_ecp_monicity_search`, substitution-map metadata, and exact stage replay for variable changes. That replay records `:not_moved` for non-first monic entries, so it cannot satisfy Park-Woodburn Section 4's first-entry precondition.

## Approach Options

1. Add an internal `ECPMonicityNormalization` record for #243 contexts and teach the existing monicity stage to replay first-coordinate elementary moves.
   This reuses the existing search, substitution, factor replay, and staged certificate patterns while adding the missing context-level API.

2. Only extend `_deterministic_ecp_monicity_search`.
   This would fix existing variable-change certificates but would not expose a context-level record, identity normalization, or staged failure diagnostics for #243 contexts.

3. Push first-entry normalization into `ecp_link_witness`.
   This would keep the link path narrow but would mix preprocessing with Park-Woodburn link witness data and still leave no reusable normalization record.

The chosen design is option 1. It is the smallest change that satisfies the new interface while preserving existing exact reduction semantics.

## Design

Add unexported structs in `src/algorithm/column_reduction.jl`:

- `ECPMonicityNormalization`, a replay-verified success record.
- `ECPMonicityNormalizationFailure`, a staged diagnostic returned when no identity or bounded variable-change candidate verifies.

Add `ecp_monicity_normalization(context::ECPInputContext; variable_order = context.variable_order, selected_variable = context.selected_variable, selected_monic_index = nothing, max_shift_power = 3, shift_signs = (one(context.ring), -one(context.ring)))`.

The constructor verifies the input context, resolves the selected variable, moves that variable to the end of the search order, and tries candidates in deterministic order:

1. Identity substitution, so already-monic #242 contexts can produce a normalization record without artificial variable changes.
2. Existing bounded shifts of one source variable by signed powers of the selected target variable.

For each candidate, the constructor:

- builds forward and inverse substitution maps;
- checks those maps are inverse polynomial-ring automorphisms on the ring generators;
- transforms the column;
- selects the first monic entry in the selected variable, or verifies a supplied `selected_monic_index` hint;
- builds elementary coordinate-move factors that send the selected monic entry to coordinate 1;
- reduces the normalized column with the existing supported exact reducer;
- inverse-substitutes the coordinate-move and reduction factors back over the original ring;
- verifies the full factor pipeline sends the original column exactly to `e_n`;
- stores replay metadata and rejects any internal mismatch.

The coordinate move for `i != 1` uses the standard elementary sequence:

```julia
E(1, i, 1), E(i, 1, -1), E(1, i, 1)
```

Applied to the transformed column, this puts entry `i` in row 1 and keeps every factor elementary over the same polynomial ring. When `i == 1`, the coordinate-move sequence is empty and the strategy is `:already_first`.

The existing `:monicity_normalization` certificate stage will use the same coordinate-move helper. Its replay will recompute the normalized column, inverse-substituted coordinate moves, inverse-substituted reduction factors, full factors, and verification tuple. A malformed stage that leaves coordinate moves empty while `selected_monic_index != 1` will fail replay.

## Tests

Add `test/expert/ecp_monicity_normalization.jl` and register it in the expert group. The focused verification command is:

```bash
julia --project=. -e 'include("test/expert/ecp_monicity_normalization.jl")'
```

The test covers:

- one #242 context that is already first-monic;
- one bounded variable-change case where a selected entry becomes monic after the search;
- one case where the selected monic entry starts outside position 1 and yields non-empty coordinate-move factors;
- inverse automorphism checks for forward and inverse substitution maps on every ring generator;
- replay of coordinate-move and inverse-substituted factors over the original ring;
- exact full-pipeline reduction of the original column to `e_n`;
- negative controls for corrupt inverse substitution, selected variable, selected monic entry, coordinate-move factors, and the empty-move-for-non-first case.

## Out Of Scope

This change does not extract Park-Woodburn link witnesses, realize `SL_3` link obligations, route through `reduce_unimodular_column`, or export new public API names.

## Self-Review

- No placeholders remain.
- The selected approach keeps the new API internal and follows existing replay-verifier patterns.
- The coordinate move is concrete elementary data, not metadata.
- Search failure returns a diagnostic record instead of guessing.
- The test cases map directly to the hard acceptance requirements.
