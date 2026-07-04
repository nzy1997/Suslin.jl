# Issue 292 Conservative Steinberg Commutator Rewrites Design

Issue #292 adds conservative internal rewrites for the Park-Woodburn Section 6
four-factor Steinberg commutator windows. The pass must consume only exact
canonical elementary-factor records, record every accepted rewrite in a
Steinberg optimization certificate, and preserve exact products.

## Context

There is no repository `AGENTS.md`; the README, existing Superpowers docs, and
the merged Steinberg optimizer issues are the working context. Issue #292 has
no GitHub comments. Dependency #291 is closed by merged PR #308, which added
the private adjacent rewrite optimizer in `src/algorithm/redundancy.jl`.
Dependency #290 is closed by PR #306 and provides
`SteinbergOptimizationCertificate`, `_steinberg_optimization_certificate`, and
`_verify_steinberg_optimization_certificate`.

The current branch already contains the #288 fixture catalog. That catalog has
positive examples for:

- `steinberg-commutator-forward-qq`;
- `steinberg-commutator-reverse-qq`;
- `steinberg-disjoint-commutator-identity-qq`.

The run is non-interactive under Agent Desk, so the Standing Answer Policy is
used for approvals. No visual companion is needed because the design has no UI
or diagram-dependent decision.

## Approach Options

Recommended: add one private commutator optimizer beside the adjacent optimizer
in `src/algorithm/redundancy.jl`. It scans canonical records left to right,
accepts only exact four-record windows, locally verifies the proposed window
rewrite by exact matrix multiplication, emits the shortened replacement, and
returns a `SteinbergOptimizationCertificate`. This keeps #292 scoped to the
new Section 6 rules without changing public APIs or default factorization
routes.

Alternative: fold commutator matching into
`_steinberg_adjacent_rewrite_optimization_certificate`. That would make one
internal pass, but it mixes two different matching granularities and makes
negative controls harder to reason about.

Alternative: add a public optimizer now. That belongs to #293 and would
prematurely freeze user-facing API behavior.

## Chosen Design

Add `_steinberg_commutator_rewrite_optimization_certificate(factors)` as a
non-exported helper. It first builds `_steinberg_sequence_context(factors,
"original")`, reusing ordinary-polynomial ring validation and canonical
factor records. It then scans the record sequence.

At each index, the pass considers only the next four records. A candidate
window is accepted only when all four records are `kind = :elementary`, share
the common ring and matrix size already validated by the sequence context, and
match one of these exact shapes:

- forward commutator:
  `E_ij(a) E_jl(b) E_ij(-a) E_jl(-b) -> E_il(a*b)`, with `i != l`;
- reverse commutator:
  `E_ij(a) E_li(b) E_ij(-a) E_li(-b) -> E_lj(-(a*b))`, with `j != l`;
- disjoint commutator identity:
  `E_ij(a) E_lp(b) E_ij(-a) E_lp(-b) -> ()`, with `i != p` and `j != l`.

The pass tries the two nontrivial commutators before the disjoint identity.
Those predicates are mutually exclusive under the stated inequalities, but the
order keeps the implementation aligned with the issue text.

Before accepting a candidate, the pass multiplies the original four matrices
and the proposed replacement matrices over the common ring and requires exact
equality. If the local product check fails, the window is skipped and the
current original factor is copied unchanged. This makes the matcher
conservative even if future coefficient or record behavior changes.

When a rewrite is accepted, the pass appends the replacement factors, records
a rewrite log entry using the existing span/count shape, and advances by four
records. When no rewrite is accepted, it appends the original factor object at
the current index and advances by one. This preserves nonmatching windows
without unrelated canonical re-materialization.

Rewrite records use rule names:

- `:commutator_forward`;
- `:commutator_reverse`;
- `:disjoint_commutator_identity`.

Metadata records the matched indices, coefficients `a` and `b`, and
`local_products_equal = true`. The final certificate is created through
`_steinberg_optimization_certificate`, so global exact product equality,
factor-count delta, metrics, and rewrite-log counts are verified by the
existing #290 machinery.

## Tests

Extend `test/expert/steinberg_factor_count_optimization.jl` with focused
tests for the new private helper:

- a hand-built forward window `E_12(a), E_23(b), E_12(-a), E_23(-b)` rewrites
  to one `E_13(a*b)` factor;
- a hand-built reverse window rewrites to one factor with coefficient
  `-(a*b)`;
- a hand-built disjoint commutator rewrites to no factors;
- each positive certificate verifies, records the expected rule name, reduces
  factor count by the expected amount, and has exact product equality;
- fixture-backed checks compare optimized factors with the #288 expected
  rewrite factors for the three commutator catalog cases;
- negative controls cover reordered factors, wrong inverse coefficients, and
  invalid index inequalities, and require the candidate sequence to stay
  unchanged with an empty commutator rewrite log.

No public API test is added because the helper remains internal and #293 owns
the public optimizer.

## Verification

Focused verification:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Full verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Both commands must exit 0. The focused command must demonstrate exact product
equality for all accepted commutator rewrites and unchanged behavior for the
negative controls.

## Out Of Scope

Do not search globally for an optimal factor count. Do not add performance
claims. Do not export a public optimizer. Do not change
`elementary_factorization(A)` or optimize any public route by default. Do not
broaden the optimizer beyond ordinary polynomial elementary factor sequences.

## Self-Review

This design has no placeholders, keeps the behavior internal, and ties every
issue requirement to either exact canonical-record matching, local product
checking, certificate verification, or focused tests. The runtime pass is
intentionally narrow: nonmatching or ambiguous windows are copied unchanged,
and public API decisions remain deferred to #293.
