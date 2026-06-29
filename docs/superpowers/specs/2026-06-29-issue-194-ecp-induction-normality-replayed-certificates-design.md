# Issue 194 ECP Induction Normality Replayed Certificates Design

## Context

Issue 194 updates the staged elementary column property induction/normality
route introduced for issue 90. The current route accepts fixture-shaped
`normality_witness` data, calls `realize_conjugate_elementary`, stores raw
rewrite factors in a named tuple, and verifies by rebuilding that same raw
rewrite. That preserves staged behavior, but it does not consume the replayable
conjugated-elementary certificate API added by issue 193.

Issue 193 is present on `main` through merge commit `d0cc533`, which added
`ConjugatedElementaryNormalityCertificate`,
`realize_conjugate_elementary_certificate`, and
`verify_conjugate_elementary_certificate`. Live GitHub issue reads are blocked
in this sandbox by the configured proxy, so the issue body supplied by Agent
Desk, the checked-in issue 90 and issue 193 specs, and the local merged #193
commits are the binding context.

## Chosen Approach

Extend `ECPInductionNormalityCertificate` with a nested normality certificate
field while keeping the existing `normality_witness` field for compatibility.
The constructor accepts either:

- the existing witness shape with `source = :supplied_normality_witness`,
  `conjugator`, `sl2_indices`, and `sl2_entry`; or
- a supplied `ConjugatedElementaryNormalityCertificate` through a new
  `normality_certificate` keyword.

When only the old witness is supplied, the route constructs the nested
certificate with
`realize_conjugate_elementary_certificate(conjugator, fixed_index, moving_index, entry)`.
The rewrite named tuple stores the nested certificate, its factors, its
product, and the existing replay metadata. The final ECP factors remain:

```text
lifted lower-variable reduction factors
nested normality certificate factors
link-step reduction factors
```

The verifier must not trust stored factor vectors. It replays the lower
reduction, rebuilds or reverifies the nested normality certificate, checks that
the stored nested certificate matches the witness data and #193 verifier, then
recomputes final factors and the final reduced column.

## Alternatives Considered

1. **Compatibility conversion plus optional supplied certificate.** Selected.
   This satisfies the new #193 replay semantics without breaking tests and
   callers that still pass `normality_witness` directly.
2. **Require callers to pass only `normality_certificate`.** Rejected because
   issue 194 explicitly asks to keep backwards compatibility for existing
   tests that pass `normality_witness`.
3. **Keep storing only the raw rewrite factors and call the #193 helper only in
   the constructor.** Rejected because the verifier would still trust stored
   raw data instead of recording and checking the nested certificate as part of
   ECP replay.

## Verification Rules

The constructor rejects a supplied nested certificate unless:

- it verifies with `verify_conjugate_elementary_certificate`;
- its `A`, `i`, `j`, and `a` match the normalized witness data;
- its ring and dimension match the ECP ring and column length;
- its convention is `:A_E_invA`; and
- its target product fixes the lower-variable column.

The ECP verifier returns `false` unless a fresh replay confirms:

- the link-step certificate still verifies;
- the lower-variable reduction still verifies or, for raw factor sequences,
  still exactly reduces `v(0)` to `e_n`;
- lifted lower factors match the replayed lower factors;
- the nested #193 certificate verifies and matches either the supplied witness
  data or the stored certificate input;
- the stored `normality_rewrite.normality_certificate` equals the replayed
  nested certificate;
- stored rewrite factors equal the nested certificate factors;
- stored final factors equal lifted lower factors plus nested factors plus
  link-step reduction factors; and
- applying the final factors to the original column gives `e_n`.

## Files

- Modify `src/algorithm/column_reduction.jl` to add the nested certificate
  field, constructor keyword, witness-to-certificate conversion, supplied
  certificate validation, and verifier replay checks.
- Modify `test/expert/ecp_induction_normality.jl` to cover constructed and
  supplied nested certificates, plus tampering of nested certificate fields.
- Modify `test/expert/elementary_column_property.jl` to assert the public
  staged route records a #193-style nested certificate.

## Out Of Scope

- No general polynomial ECP reducer from issue 185.
- No Murthy, Quillen, recursive `SL_n`, Laurent/ToricBuilder, or
  Steinberg factor-count optimization.
- No change to the public `reduce_unimodular_column` API.

## Verification

Focused commands:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused tests must include an ECP induction/normality certificate whose
nested normality field is a `ConjugatedElementaryNormalityCertificate`, whose
final factors still reduce the staged column exactly as before, and whose
verifier rejects tampering with a child factor or stored target matrix inside
the nested proof.

## Spec Self-Review

- The design preserves existing witness compatibility.
- The nested replay object is the #193 conjugated-elementary certificate, not a
  new parallel normality format.
- Verification recomputes or reverifies the nested certificate and final factor
  sequence.
- The scope is limited to the staged ECP induction/normality path.
