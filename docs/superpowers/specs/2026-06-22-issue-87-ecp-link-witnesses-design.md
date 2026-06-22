# Issue 87 ECP Link Witnesses Design

## Context

Issue #87 records the Park-Woodburn Elementary Column Property link-theorem
data that comes after a column has a selected monic first entry. Existing ECP
work already provides fixture-backed ordinary polynomial columns, replayable
column reduction certificates, variable-change replay metadata, and deterministic
monicity search records. This issue must not implement the link-step operations
or return reducer factors.

## Chosen Approach

Add a thin, non-exported expert witness record in `src/algorithm/column_reduction.jl`
and focused expert tests in `test/expert/ecp_link_witnesses.jl`.

The record accepts supplied link-theorem data for ordinary polynomial columns
whose first entry is monic in the selected variable. The verifier recomputes
every recorded algebraic identity:

- the selected first entry equals `v[1]` and is monic in the selected variable;
- each tail reduction has lifted tail coefficients whose exact linear
  combination of `v[2:end]` gives the recorded `tilde_G`;
- each recorded resultant equals `resultant(v[1], tilde_G, selected_variable)`;
- each Bezout pair satisfies `f*v[1] + h*tilde_G == resultant`;
- the coverage multipliers satisfy `sum(resultant_i*g_i) == 1`;
- the path starts at `0`, ends at the selected variable, and advances by
  `resultant_i*g_i*selected_variable`.

Extraction is deliberately staged. Calling the construction helper without
`supplied_link_witness` throws `ArgumentError` that names the missing supplied
or extracted Park-Woodburn link witness data. Supplied records carry
`:supplied_link_witness` metadata so tests mark the boundary explicitly.

## Alternatives Considered

1. **Supplied witness replay first.** This is the chosen route. It satisfies the
   issue's validation objective, reuses existing fixture IDs, and keeps the
   later operation stage separate.
2. **Implement extraction with Oscar Groebner/resultant helpers now.** Deferred:
   the issue permits supplied witnesses, and extraction would broaden the scope
   beyond fixture-backed validation.
3. **Attach link data to `ECPColumnReductionCertificate`.** Deferred: #84's
   guardrail asks for a thin replay layer and no broad public schema. A separate
   non-exported witness record avoids changing public reducer behavior.

## Interfaces

Production code adds:

- `ECPLinkWitnessRecord`
- `ecp_link_witness(v, R; variable_order, selected_variable, selected_monic_index,
  supplied_link_witness)`
- `verify_ecp_link_witness(record)::Bool`

These names remain non-exported, matching the existing expert/internal ECP APIs.
Expert tests access them as `Suslin.<name>`.

The supplied witness is a `NamedTuple` with exact metadata:

- `source = :supplied_link_witness`
- `residue_probes`
- `tail_reductions`
- `resultants`
- `bezout_coefficients`
- `coverage_multipliers`
- `path_points`

The constructor validates immediately and throws `ArgumentError` on corrupted
or incomplete supplied data. The verifier returns `false` for tampered stored
records.

## Fixtures And Tests

The expert test reuses #83 fixture IDs:

- `ecp-variable-change-monic-gf2`
- `ecp-monic-first-entry-qq`

Both are ordinary polynomial columns with no unit entry and a monic first entry
in the selected variable. Tests supply exact Park-Woodburn link witness data for
each fixture, including one one-resultant GF(2) cover and one two-resultant QQ
cover.

Negative controls corrupt one resultant, one Bezout coefficient, one coverage
multiplier, and one path point. Each corrupted construction throws, and tampered
records fail `verify_ecp_link_witness(record)` without returning reducer
factors.

## Out Of Scope

- No link-step operation from `v(b_{i-1})` to `v(b_i)`.
- No public `reduce_unimodular_column` return-type change.
- No Quillen patching or final Park-Woodburn matrix driver.
- No claim that ordinary witness-unit data is enough for the link theorem.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witnesses.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

The existing expert harness should also register the new expert file.
