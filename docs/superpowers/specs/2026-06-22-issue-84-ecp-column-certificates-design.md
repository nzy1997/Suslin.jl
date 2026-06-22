# Issue 84 ECP Column Certificates Design

## Context

`reduce_unimodular_column(v, R)` reduces currently supported unimodular
columns to `e_n` and returns only the elementary factor sequence. Issue 84
needs an expert replay path for those existing reducer families so later
Elementary Column Property work can inspect validation, variable changes,
monicity forcing, ideal witnesses, embedded operations, and Laurent
normalization without changing the public reducer return type.

Issue 83 is merged and provides an ECP column fixture catalog covering direct
unit, witness-unit, variable-change monicity, embedded block, staged failure,
and link-oriented metadata examples. Issue 70 and PR 77 provide the local
`SL_3` certificate precedent: internal certificate names are acceptable, but
every recorded field must participate in exact replay.

The issue comment is binding: keep the certificate layer thin and
replay-driven. Stage tags may name the route, but verification must not accept a
certificate because a tag is present.

## Design Choice

Add a non-exported expert certificate path inside
`src/algorithm/column_reduction.jl`:

- `ECPColumnReductionCertificate`
- `ecp_column_reduction_certificate(v, R)`
- `verify_ecp_column_reduction(cert)::Bool`

The legacy `reduce_unimodular_column(v, R)` delegates to the same certificate
construction path and returns `certificate.factors`. The new names remain
available as `Suslin.<name>` for expert tests, but they are not exported from
`src/Suslin.jl`.

Alternatives considered:

- Change `reduce_unimodular_column` to return a rich result. This violates the
  issue objective and would break existing callers.
- Add a parallel reducer just for certificates. This duplicates stage logic and
  risks drift between factor-returning and certificate-returning behavior.
- Record broad fixture metadata from issue 83 wholesale. That would create
  decorative fields. The certificate records only data the verifier checks.

## Certificate Shape

The certificate records:

- `original_column`: the coerced input column over `R`.
- `ring`: the target coefficient ring.
- `stages`: a tuple of replayable stage records.
- `factors`: the final factor sequence over the original ring.
- `final_column`: the stored result of applying `factors` to
  `original_column`.
- `verification`: the verification summary computed at construction time.

Current stage kinds are:

- `:validation`: records the input length and unimodularity flag. Replay checks
  the column length and recomputes `is_unimodular_column`.
- `:unit_entry`: records the pivot index and pivot unit. Replay checks the
  pivot value is unchanged, unit/inverse relation holds, recomputes the exact
  unit-entry factors, and checks their product sends the stage input to `e_n`.
- `:witness_unit`: records the exact ideal-membership witness, pivot index, and
  witness unit. Replay checks `sum(w[i] * v[i]) == 1`, the witness pivot is a
  unit, recomputes unit-creation and unit-reduction factors, and checks the
  stage output.
- `:monicity_normalization`: records the selected variable index, last
  variable, shift power/sign, forward and inverse substitution values,
  substituted column, transformed factors, and inverse-substituted factors.
  Replay checks substitution exactly, checks a transformed entry is monic in the
  last variable, replays the transformed unit/witness stage, checks inverse
  substitution of each factor, and checks the inverse-substituted factors reduce
  the original column.
- `:embedded_three_block`: records selected indices, subcolumn, subcertificate,
  embedded block factors, post-block column, elimination factors, and move
  factors. Replay checks the subcertificate, block embedding, post-block column,
  elimination coefficients, move factors, and final column.
- `:laurent_normalization`: records the Laurent normalization object, normalized
  polynomial column, polynomial subcertificate, lifted factors, normalization
  unit factors, and shift units. Replay checks normalization metadata, the
  polynomial subcertificate, lifted factors, unit normalization factors, and the
  final Laurent reduction.

Stages can be nested only where the existing reducer is nested: monicity
contains the transformed unit/witness stage, embedded block contains the
three-entry subcertificate, and Laurent normalization contains the polynomial
subcertificate. Unsupported and non-unimodular cases continue to throw the
existing errors and do not produce certificates.

## Replay Rules

`verify_ecp_column_reduction(cert)` returns `false` for malformed or tampered
certificates, except `InterruptException` is rethrown. Verification:

1. Checks the certificate has a valid ring, column, factor vector, and target
   shape.
2. Replays each stage from its recorded input, never trusting the stage tag
   alone.
3. Recomputes witness equations, substitutions, embedded factors, lifted
   factors, and stage products from the recorded algebraic data.
4. Checks the concatenated replay factor sequence equals `cert.factors`
   element by element.
5. Checks `cert.final_column` equals `product(cert.factors) *
   original_column`.
6. Checks the final column is exactly `e_n`.

The verifier rejects changed factors, changed selected variables or inverse
substitutions, changed witness coefficients, changed embedded indices, changed
Laurent shift units, and any extra or missing witness keys for stage records
that have a fixed witness shape.

## Tests

Add `test/expert/ecp_column_certificate.jl` and register it in the expert
group. The focused test builds certificates for the current successful reducer
families:

- direct unit entry from `ecp-unit-entry-gf2`,
- witness-unit from `ecp-witness-unit-gf2`,
- monicity normalization from `ecp-variable-change-monic-gf2`,
- embedded three-entry block from `ecp-longer-embedded-block-gf2`,
- Laurent-normalized existing coverage from
  `test/expert/unimodular_reduction_exact.jl`.

For each certificate the test checks:

- `verify_ecp_column_reduction(cert) == true`,
- `cert.factors` reduce `cert.original_column` to `e_n`,
- `reduce_unimodular_column(v, R)` still returns a factor vector reducing the
  same original column exactly.

Negative controls mutate one recorded factor, one monicity inverse
substitution, one witness coefficient, and one embedded-block index or Laurent
shift field. Each tampered certificate must return `false` or throw a staged
verification error. The test must also prove unsupported and non-unimodular
fixture cases keep their existing staged errors.

## Verification

Focused expert command:

```bash
julia --project=. -e 'include("test/expert/ecp_column_certificate.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Full expert-inclusive command:

```bash
julia --project=. test/runtests.jl expert
```

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to the existing successful reducer families.
- The public reducer still returns only factors.
- Each certificate field has a stated replay check.
- Unsupported columns remain unsupported.
