# Issue 195 Polynomial Normality Support Boundary Design

## Context

Issues 190 through 194 are merged on `main` through commit `708131d`, so the
ordinary-polynomial normality layer now has fixture-backed tests for:

- Cohn-type replay certificates (`test/expert/cohn_type.jl`);
- orthogonal rank-one replay certificates (`test/expert/normality_rank_one.jl`);
- conjugated-elementary replay certificates (`test/expert/normality.jl`); and
- the ECP induction adapter that stores a nested conjugated-elementary
  certificate (`test/expert/ecp_induction_normality.jl`).

Issue 195 is a closeout gate for parent issue 181. It must make that support
boundary visible in public docs and tests without implying support for the
later mainline algorithm issues 182 through 187 or issue 188 optimization work.

GitHub REST issue reads were available for issue 195 and the merged dependent
pull requests. GitHub GraphQL and browser reads were unavailable in this
sandbox because the configured proxy denied the connection, so the issue body,
issue comment, REST pull request bodies, and checked-in code are the binding
context.

## Chosen Approach

Add a narrow documentation boundary plus an expert support-boundary smoke gate.
The README and Documenter index will name the supported ordinary-polynomial
certificate layer:

- Cohn-type realization certificates;
- rank-one normality certificates;
- conjugated-elementary normality certificates; and
- the staged ECP induction/normality adapter that replays the nested
  conjugated-elementary certificate.

The same docs must state that arbitrary Park-Woodburn `SL_n(k[x_1, ..., x_m])`
factorization, Murthy local solving, full Quillen patching, the general ECP
reducer, recursive `SL_n` mainline acceptance, Laurent/ToricBuilder acceptance,
and Steinberg factor-count optimization remain staged work.

The test gate will live in a focused expert file:
`test/expert/polynomial_normality_support_boundary.jl`. It will scan
`README.md` and `docs/src/index.md` for the supported ordinary-polynomial
normality certificate wording and for explicit staged-boundary wording. It will
also assert that the expert test registry includes the four certificate suites
introduced by issues 190 through 194. Existing public negative controls remain
the behavioral guard for unsupported full Park-Woodburn claims.

## Alternatives Considered

1. **Docs plus explicit expert support-boundary gate.** Selected because it
   makes the public boundary testable and keeps issue 195 out of algorithm
   implementation.
2. **Only update README and docs.** Rejected because the issue asks for a test
   gate proving the normality certificate suite remains visible.
3. **Add new public API wrappers or algorithm acceptance.** Rejected because
   that would cross into issues 182 through 187 and overstate support.

## Validation Rules

The expert support-boundary gate must verify:

- both README and docs index name ordinary-polynomial normality/conjugation
  certificates;
- both docs name Cohn-type, rank-one, and conjugated-elementary certificates;
- both docs mention the staged ECP nested-normality adapter;
- both docs preserve a staged boundary for full Park-Woodburn
  `SL_n(k[x_1, ..., x_m])` and later algorithm work; and
- `test/runtests.jl` expert registration includes `cohn_type.jl`,
  `normality_rank_one.jl`, `normality.jl`, and
  `ecp_induction_normality.jl`.

The implementation must not change algorithm behavior, public factorization
acceptance, Laurent/ToricBuilder routes, or certificate APIs.

## Files

- Modify `README.md`.
- Modify `docs/src/index.md`.
- Modify `test/runtests.jl`.
- Add `test/expert/polynomial_normality_support_boundary.jl`.

## Verification

Focused red/green command:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_support_boundary.jl")'
```

Required issue commands:

```bash
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
```

The expert suite must include the Cohn-type, rank-one, conjugated-elementary,
and ECP nested-normality checks. The default package entry point must remain the
fast `public` and `internal` suites.

## Spec Self-Review

- No placeholders remain.
- The scope is documentation and test registration only.
- The support wording is limited to ordinary-polynomial normality/conjugation
  certificates.
- The staged-failure boundary remains explicit for full Park-Woodburn,
  Laurent/ToricBuilder, and Steinberg optimization work.
