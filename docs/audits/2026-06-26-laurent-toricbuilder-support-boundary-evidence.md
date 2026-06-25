# Laurent/ToricBuilder Support-Boundary Evidence

Date: 2026-06-26

This page reconciles the support-scope audit in #129 with the ToricBuilder
Q-block status in #131 after the `case_010` Laurent certificate route (#135)
and the bounded `case_008` triage (#137). It is evidence for review only; it
does not close the parent issues and does not add algorithm support.

## Support Matrix

| Evidence item | Route | Supported outcome | Boundary |
| --- | --- | --- | --- |
| `case_010` ToricBuilder Q-block | `laurent_gl_factorization_certificate` | `gl_certificate_pass`; verified `true`; decomposed base matrices `48` | public `elementary_factorization` remains `staged_boundary` for the original Laurent `GL_n` input |
| `case_008` bounded exercise | bounded Laurent `GL_n` certificate route | `certified_algorithm_boundary` at `certificate_construction` under explicit `--exercise=case_008 --timeout-seconds=120` | not a default report pass; remains a staged algorithm boundary |
| Default ToricBuilder Q-block report rows `case_001`-`case_006`, `case_010` | generated Q-block status report | `gl_certificate_pass` with public `staged_boundary` and Laurent monomial-unit determinants | evidence is the Laurent `GL_n` certificate route, not original-input elementary factor sequences |
| Laurent monomial-unit `GL_n` inputs in the staged certificate path | `normalize_laurent_gl_matrix` then `laurent_gl_factorization_certificate` | certificate verifies the normalized determinant-one core and exact reconstruction metadata | original-input `elementary_factorization` for Laurent `GL_n` remains a `staged boundary` |
| Ordinary-polynomial staged slices | `elementary_factorization` | exact elementary factor sequences that satisfy `verify_factorization(A, factors) == true` | not arbitrary Park-Woodburn support |
| Remaining Laurent scope | no broad public factor-sequence route | certificate-backed monomial-unit slices only where recorded tests exercise them | not arbitrary Laurent `GL_n` support |

## Route Diagram

`elementary_factorization(A) -> exact elementary factor sequence -> verify_factorization(A, factors)`

This is the public factor-sequence route for the supported ordinary-polynomial
and determinant-one Laurent `SL` staged slices. It returns factors only when the
current staged implementation can verify exact multiplication back to `A`.

`ToricBuilder Q-block -> classify Laurent determinant -> normalize Laurent GL_n determinant -> factor determinant-one core -> verify Laurent GL_n certificate`

This is the Laurent `GL_n` certificate route. It records determinant
normalization, factors the normalized determinant-one core, and verifies the
certificate metadata. It is not the same as returning an elementary factor
sequence for the original Laurent `GL_n` input.

## Verification Commands

Run the issue-required documentation and algorithm checks:

```bash
julia --project=. test/runtests.jl
julia --project=. test/runtests.jl expert
```

Run the package entry point used by Agent Desk workers:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Optional route evidence commands:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_010 --output=/tmp/case010-q-block-status.md
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=120 --output=/tmp/qblock-case008.md
```

Expected `case_010` evidence: the generated row contains
`gl_certificate_pass`, `verified true`, public `staged_boundary`, and a positive
decomposed base-matrix count.

Expected `case_008` evidence: the bounded generated row is structured, not
`not_exercised_in_default_report`, not `route_error`, and not an unstructured
timeout. The current observed bounded outcome is `certified_algorithm_boundary`
at `certificate_construction`.

## Parent-Issue Reconciliation

- #129 remains the support-scope audit. This page confirms the staged support
  boundary: exact factor sequences are available only on the supported
  `elementary_factorization` slices, while Laurent `GL_n` monomial-unit support
  is certificate evidence unless a later issue adds original-input factor
  sequences.
- #131 remains the ToricBuilder Q-block status thread. The default generated
  report now records `case_010` as a Laurent `GL_n` certificate pass, while
  `case_008` is documented through the explicit bounded exercise route.
- #135 is reflected as the supported `case_010` certificate outcome.
- #137 is reflected as the bounded `case_008` structured boundary outcome.

## Explicit Non-Claims

- This is not arbitrary Park-Woodburn `SL_n(k[x_1, ..., x_m])` support.
- This is not arbitrary Laurent `GL_n` support.
- This is not a claim that `elementary_factorization` returns factor sequences
  for original Laurent `GL_n` inputs with monomial-unit determinant.
- This is not a performance claim or benchmark table.
