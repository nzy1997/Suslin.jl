# Issue 86 ECP Monicity Search Design

## Context

Issue #86 replaces the current hard-coded monicity normalization loop with a deterministic helper that is explicit about the selected variable order, searched shifts, accepted variable change, and staged failure when the bounded search is exhausted. The implementation builds on #83 fixture catalog entries, #84 replayable ECP reduction certificates, and #85 variable-change stage records.

The repository has no `AGENTS.md`. The current worker branch starts from `main` at the merge of PR #93, so the #85 record fields are available in `src/algorithm/column_reduction.jl`.

## Chosen Approach

Use a thin internal search-result layer around the existing #85 variable-change stage builder.

This approach keeps `reduce_unimodular_column(v, R)` returning factors, preserves the certificate record shape for successful monicity normalization, and adds a replayable expert-facing result or failure object for the deterministic search itself. The helper searches in a supplied variable order, chooses a target last variable from that order, then tries earlier source variables and bounded signed powers in stable nested-loop order. On success it returns the existing monicity stage data plus the accepted search metadata. On exhaustion it returns a staged failure naming the variable order, source variables, target variable, shift powers, signs, and reducer tried.

## Alternatives Considered

1. Replace the reducer stage directly with a public monicity-search API. This would expose more surface than #86 needs and would conflict with the current pattern where expert tests access non-exported internals through `Suslin.<name>`.
2. Keep the old loop and only add failure messages. This would not satisfy the variable-order requirement or provide a replayable search record.
3. Add a thin helper and route the existing stage through it. This is selected because it satisfies #86 while reusing #85 replay fields and keeping unsupported cases staged.

## Interface

Add internal record structs in `src/algorithm/column_reduction.jl`:

- `ECPMonicitySearchResult`: contains the original column, ring, requested variable order, bound, target/source variable metadata, shift power/sign/polynomial, transformed column, selected monic entry/index, and the successful #85 monicity stage.
- `ECPMonicitySearchFailure`: contains the original column, ring, requested variable order, bound, search space, attempted candidate count, and a message that explicitly says the bounded deterministic monicity search was exhausted.

Add expert helper:

- `_deterministic_ecp_monicity_search(column, R; variable_order = tuple(gens(R)...), max_shift_power = 3, shift_signs = (one(R), -one(R)))`

The helper returns either `ECPMonicitySearchResult` or `ECPMonicitySearchFailure`. It must not throw for a unimodular-but-exhausted ordinary polynomial column, and it must not claim the column is non-unimodular.

## Reducer Integration

`_reduce_after_monicity_normalization_certificate` will call the helper with the ring generator order and default bound. On success it returns the helper's staged #85 monicity record and factors. On failure it returns `nothing` so the existing higher-level unsupported error remains the public behavior.

## Failure Semantics

The staged failure object records:

- `kind = :monicity_search_exhausted`
- the requested variable order as actual ring generators
- `max_shift_power`
- target variable and source variables
- shift signs and candidate shift polynomials
- the number of attempted candidates
- a stable message suitable for tests

A valid but incomplete variable order is part of the bounded search space. If the supplied order omits the source or target variable needed by a supported fixture, the helper returns `ECPMonicitySearchFailure`. Invalid search inputs such as duplicated variables, non-generator variables, or a negative bound throw `ArgumentError` because they are caller errors rather than staged mathematical failure. A zero bound is allowed and means no positive-power shifts are searched.

## Testing

Add `test/expert/ecp_monicity_search.jl` with fixture-backed tests:

- The existing `ecp-variable-change-monic-gf2` case succeeds with the default variable order and produces a replayable result whose stage reduces exactly to `e_n`.
- The existing `ecp-variable-change-permuted-gf2` case succeeds only when the deterministic order allows the required target/source choice, proving the search is driven by the supplied order rather than only the ring's native generator order.
- The `ecp-unsupported-unimodular-gf2` case staged-fails for a deliberately exhausted bounded search, with variable order and bound present in the failure.
- Negative controls lower the bound or remove the required variable and assert staged failure or caller validation instead of unrecorded substitution.

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_monicity_search.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out of Scope

This issue does not add randomness, Quillen patching, local-to-global patching, or support for arbitrary Noether-normalization choices. It does not change the public return type of `reduce_unimodular_column`.

## Self-Review

- No placeholder sections remain.
- The design is scoped to `src/algorithm/column_reduction.jl` and `test/expert/ecp_monicity_search.jl`.
- The selected helper keeps the #85 record fields replayable and preserves existing public behavior.
- Staged failure covers bounded search exhaustion, including valid orders that omit the needed variable; malformed search requests remain explicit caller errors.
