# Issue 88 ECP Link Step Design

## Context

Issue #88 builds on the supplied Park-Woodburn link-theorem witness records from
#87. Those records already replay the selected monic first entry, tail
reductions, resultants, Bezout identities, coverage multipliers, and path points
`b_0 = 0, ..., b_L = X`. This issue adds the next supported operation stage:
turn a validated link witness into a replayable certificate that relates each
path column `v(b_{i-1})` to `v(b_i)` and records the lower-variable obligation
`v(0)`.

## Chosen Approach

Add a narrow internal `ECPLinkStepCertificate` in
`src/algorithm/column_reduction.jl` and fixture-backed expert coverage in
`test/expert/ecp_link_step.jl`.

The constructor accepts a validated `ECPLinkWitnessRecord`, evaluates each path
column by substituting the selected variable with the recorded path point, and
builds one segment record per path step. Each segment replays the #87 identity
for that exact step:

- `delta = b_i - b_{i-1}`;
- `delta == selected_variable * resultant_i * coverage_multiplier_i`;
- the tail coefficients recompute `tilde_G_i`;
- the Bezout pair at `b_{i-1}` proves the step ideal identity;
- `v(b_i) - v(b_{i-1})` is exactly divisible by `delta`.

For the first supported family, the implementation is deliberately limited to
the fixture-backed supplied witnesses from #87. The supported probe signatures
are replayed explicitly before construction; other verified link witnesses stage
fail instead of falling through to a generic shortcut. For those supported
fixtures, segment transport is built from existing replayable ECP
column-reduction certificates for the two endpoint path columns. If
`F_from * v(b_{i-1}) = e_n` and `F_to * v(b_i) = e_n`, the segment transport is
`F_to^{-1} * F_from`, recorded as elementary `E_n` factors and verified by
replay. The `SL_2` contribution is recorded explicitly as the supported identity
block for this family and verified as determinant one; fuller Park-Woodburn
`SL_2` block extraction remains staged for later witness families.

## Alternatives Considered

1. **Endpoint-certificate transport with exact link-identity replay.** Selected.
   It gives a conservative supported family, produces real replayable factors,
   and rejects unsupported path columns instead of inventing unchecked factors.
2. **Require caller-supplied link-step factors.** Rejected for this issue because
   #88 asks the operation to consume the #87 witness and construct the step data,
   not accept an operation certificate without replay.
3. **Implement the full paper-level `SL_2` extraction now.** Deferred. The issue
   explicitly prefers a narrow supported family over ad hoc broad support, and
   the current fixtures can be verified exactly without routing the public
   reducer or implementing lower-variable induction.

## Interface

Production code adds non-exported expert/internal names:

- `ECPLinkStepCertificate`
- `ecp_link_step_certificate(v, R; link_witness, variable_order,
  selected_variable, selected_monic_index, supplied_link_witness)`
- `verify_ecp_link_step_certificate(certificate)::Bool`

`link_witness` may be supplied directly as a validated `ECPLinkWitnessRecord`.
If it is omitted, the constructor builds one using the existing
`ecp_link_witness` supplied-witness interface. Laurent rings, malformed witness
records, unsupported path columns, unknown supplied-witness probe signatures,
and non-identity `SL_2` family requests fail with `ArgumentError`.

The certificate stores:

- original transformed column `v(X)`;
- lower-variable column `v(0)`;
- all path points and path columns;
- per-segment link identity replay, `SL_2` block contribution, elementary
  factors, forward mapping `v(b_{i-1}) -> v(b_i)`, and inverse mapping;
- composed forward factors mapping `v(0)` to `v(X)`;
- composed reduction factors mapping `v(X)` to `v(0)`;
- verification metadata for every segment and the composed maps.

## Testing

Add `test/expert/ecp_link_step.jl` with the two #87 fixture witnesses:

- `ecp-variable-change-monic-gf2`, with one nonzero step `0 -> x`;
- `ecp-monic-first-entry-qq`, with two path segments `0 -> y^2*x -> x`.

The verifier checks every path column, every segment map, the composed forward
map, the composed inverse map to the lower-variable obligation, and stored
verification equality. Negative controls tamper one segment `SL_2` block, one
segment elementary contribution, one link witness identity, and one
verified-but-unsupported supplied witness family; construction or replay
verification must fail.

## Out Of Scope

- No public `reduce_unimodular_column` routing change.
- No lower-variable induction or normality rewrite.
- No Quillen patching, random search, or factor-count optimization.
- No claim that arbitrary link witness families have non-identity `SL_2`
  extraction support yet.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```
