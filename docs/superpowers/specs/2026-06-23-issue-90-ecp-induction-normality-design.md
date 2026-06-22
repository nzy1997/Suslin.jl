# Issue 90 ECP Induction Normality Design

## Context

Issue #88 added a replayable ECP link-step certificate. That certificate verifies
the supplied Park-Woodburn path data and records elementary factors that map the
ordinary-polynomial column `v(X)` down to the lower-variable obligation `v(0)`.
Issue #90 adds the next internal stage: consume a verified link-step certificate,
consume a verified reduction of `v(0)`, replay an explicit normality/conjugation
witness, and return elementary factors over the original ring that reduce `v(X)`
to `e_n`.

The public `reduce_unimodular_column` entry point remains unchanged. This is an
expert/internal certificate path for the supported fixture-backed families only.

## Chosen Approach

Add a non-exported `ECPInductionNormalityCertificate` in
`src/algorithm/column_reduction.jl` with a constructor
`ecp_induction_normality_certificate(v, R; link_step, lower_reduction,
normality_witness)` and a verifier
`verify_ecp_induction_normality_certificate(certificate)::Bool`.

The constructor verifies the supplied #88 `ECPLinkStepCertificate`, verifies or
constructs a lower-variable reduction for `link_step.lower_variable_column`, and
lifts those lower factors into the original ring. For the first supported
normality family, the caller supplies explicit fixture-backed witness data:

- `source = :supplied_normality_witness`;
- a conjugator matrix equal to the inverse product of the lifted lower-variable
  factors;
- an embedded two-row `SL_2` elementary contribution, recorded by indices
  `(n, j)` and a nonzero entry `a`.

The replay uses `realize_conjugate_elementary(conjugator, n, j, a)` from
`src/algorithm/normality.jl` to rewrite the conjugated `SL_2` contribution as
elementary factors. Since the conjugator maps `e_n` back to `v(0)`, the
conjugated contribution fixes `v(0)`. The final factor sequence is:

```text
lifted lower-variable reduction factors
normality rewrite factors
link-step reduction factors
```

Under the existing left-multiplication convention, this first maps `v(X)` to
`v(0)`, applies a replayed normality factorization that fixes `v(0)`, and then
applies the lifted lower-variable factors to reach `e_n`. The verifier
recomputes every part, including the normality product and the final reduction
of the original column.

## Alternatives Considered

1. **Explicit normality witness plus exact final replay.** Selected. It keeps
   fixture support narrow, uses the existing normality helper, and makes the
   theorem move visible without changing public reducer routing.
2. **Only compose lower-variable factors with link-step inverse factors.**
   Rejected because it would hide the normality step and could pass without any
   recorded non-identity `SL_2` contribution.
3. **Extend #88 link-step extraction to derive arbitrary non-identity `SL_2`
   blocks.** Deferred. The issue permits deterministic fixture-backed
   normality data; arbitrary extraction belongs to a broader Park-Woodburn
   driver.

## Interface

Production code adds these non-exported expert/internal names:

- `ECPInductionNormalityCertificate`
- `ecp_induction_normality_certificate(v, R; link_step, lower_reduction,
  normality_witness)`
- `verify_ecp_induction_normality_certificate(certificate)::Bool`

`link_step` is required and must verify with `verify_ecp_link_step_certificate`.
`lower_reduction` may be omitted, in which case the existing
`ecp_column_reduction_certificate` is used for `v(0)`, or supplied as an
`ECPColumnReductionCertificate` or a concrete factor sequence. Any supplied
factor sequence must exactly reduce `v(0)` to `e_n`.

`normality_witness` is required for this stage. Unsupported witness shapes,
identity `SL_2` entries, mismatched conjugators, failed normality factorizations,
or final replay failures throw `ArgumentError` instead of returning partial
factors.

The certificate stores:

- the replayed #88 link-step certificate;
- `lower_variable_column`;
- lower-variable factors and their lifted copies in the original ring;
- the explicit normality witness and replayed normality rewrite data;
- final elementary factors;
- final column and verification metadata proving
  `product(final_factors) * original_column == e_n`.

## Testing

Add `test/expert/ecp_induction_normality.jl` and register it in the expert test
group. The test reuses the two supported #88 fixture witnesses:

- `ecp-variable-change-monic-gf2`;
- `ecp-monic-first-entry-qq`.

For each case, the test constructs and verifies a #88 link-step certificate,
constructs and verifies a lower-variable reduction certificate, supplies an
explicit non-identity normality witness, verifies the normality rewrite, and
checks that the final factors reduce the original column to `e_n`. At least one
assertion checks that the recorded `SL_2` block is non-identity.

Negative controls rebuild the stored certificate with a tampered lifted
lower-variable factor and with a tampered normality rewrite factor. Both must
make `verify_ecp_induction_normality_certificate` return `false`.

## Out Of Scope

- No public `reduce_unimodular_column` routing change.
- No Quillen local-to-global patching.
- No full Park-Woodburn matrix driver.
- No Laurent `GL_n` determinant correction.
- No factor-count optimization.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```
