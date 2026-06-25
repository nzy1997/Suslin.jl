# Issue 138 Laurent/ToricBuilder Evidence Design

## Goal

Publish one audit-style markdown page that reconciles the support-scope audit
from #129 and the ToricBuilder Q-block status from #131 after the `case_010`
certificate route and `case_008` bounded classification work.

## Context

Issue #129 records the repository's current support boundary: Suslin.jl supports
staged, certificate-backed slices of exact ordinary polynomial and Laurent
polynomial factorization, not arbitrary Park-Woodburn `SL_n(k[x_1, ..., x_m])`
factorization and not arbitrary Laurent `GL_n` elementary factorization.

Issue #131 records the ToricBuilder cache Q-block status. Since then, #135
promoted `case_010` to a verified original-input Laurent `GL_n` certificate
pass while preserving `elementary_factorization` as a public staged boundary for
that original Laurent `GL_n` input. Issue #137 classified explicit bounded
`case_008` runs as a structured bounded route outcome instead of leaving them as
an unstructured slow-case probe.

The existing generated Q-block status report remains useful as raw route data.
This issue needs a concise delivery evidence page that explains the boundary in
terms a reviewer can audit without reading the whole roadmap.

## Approach Options

Recommended: add a single hand-written evidence page under `docs/audits/` and
extend the existing expert documentation smoke test so it pins the page's
required claims. This keeps the delivery note stable and factual while the
generated Q-block status report remains the source of per-run timing rows.

Alternative: extend the generated Q-block report to include all parent-issue
reconciliation text. This would mix generated timing data with editorial audit
scope and make the report script responsible for roadmap reconciliation.

Alternative: update only README scope text. This would make the top-level README
too dense and would not give #129 and #131 a focused delivery artifact.

## Chosen Design

Create `docs/audits/2026-06-26-laurent-toricbuilder-support-boundary-evidence.md`.
The page will contain:

- a short scope statement that says this is evidence for #129 and #131 after
  #135 and #137;
- a support matrix with rows for ordinary polynomial staged factor sequences,
  ToricBuilder default Laurent Q-block certificate rows including `case_010`,
  explicit bounded `case_008`, Laurent `GL_n` certificate support, and the
  original-input `elementary_factorization` staged boundary;
- a text route diagram showing `elementary_factorization` factor sequences as a
  separate public route from Laurent `GL_n` certificate construction;
- exact verification commands from the issue plus the package test command;
- parent-issue reconciliation notes that leave #129 and #131 open for review
  and do not claim arbitrary Park-Woodburn or arbitrary Laurent `GL_n` support.

The page will use "staged boundary" exactly for the remaining original-input
`elementary_factorization` boundary. It will state that `case_010` is supported
through `laurent_gl_factorization_certificate`, not through public
`elementary_factorization` factors for the original Laurent `GL_n` input. It
will state that `case_008` is a bounded structured outcome, not a full default
certificate pass.

## Tests

Extend `test/expert/documentation_smoke.jl` with a documentation smoke test that
reads the new evidence page and checks:

- the evidence page exists;
- the support matrix names the supported `case_010` outcome;
- the support matrix names the bounded `case_008` outcome;
- the page names the Laurent `GL_n` certificate route;
- the page names the remaining staged boundary for original-input
  `elementary_factorization`;
- the page includes the exact issue verification commands;
- the page includes the conservative "not arbitrary Park-Woodburn" and "not
  arbitrary Laurent `GL_n`" scope statements.

The negative controls from the issue are covered because removing the required
support-matrix rows or replacing the exact original-input
`elementary_factorization` `staged boundary` row with a broader support claim
will remove strings the smoke test requires.

## Verification

Required issue verification:

```bash
julia --project=. test/runtests.jl
julia --project=. test/runtests.jl expert
```

Required package verification for this Agent Desk run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not add new algorithm support. Do not open or close #129 or #131. Do not
claim arbitrary Park-Woodburn support, arbitrary Laurent `GL_n` support, or
original-input Laurent `GL_n` elementary factor sequences.

## Automatic Decisions

- Clarifying questions and design approvals were answered by the Standing
  Answer Policy because this is a non-interactive Agent Desk run.
- The visual companion was skipped because this is an audit markdown/test change
  with no visual design choice that would be clearer in a browser.
- The hand-written evidence page approach was selected because it keeps support
  reconciliation separate from the generated timing report.
- The design was approved automatically under the Standing Answer Policy because
  it is the conservative approach recommended by the issue text.
