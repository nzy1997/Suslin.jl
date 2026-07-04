# Issue 290 Steinberg Optimization Certificate Design

Issue #290 adds an internal certificate and verifier for Steinberg factor
sequence rewrites. This is a replay and audit layer only: it records original
and candidate optimized elementary factor sequences, applied rewrite metadata,
before/after factor counts and PR #179 metric values, exact products, and a
verification status. It does not shorten factors or export a public API.

## Context

There is no repository `AGENTS.md`; the README, existing Superpowers docs, and
the existing private certificate patterns apply. The working branch is already
an isolated Agent Desk worktree. Issue #289 is closed by PR #304 and the branch
contains the private `_canonical_elementary_factor_record` and
`_elementary_factor_record_matrix` helpers. Issue #288 is closed by PR #302 and
the branch contains `SteinbergOptimizationFixtureCatalog` for ordinary
polynomial positive and negative Steinberg rewrite examples. PR #179 added the
metric helpers `max_elementary_factor_monomial_degree(factors)` and
`total_elementary_factor_offdiagonal_monomials(factors)`, which this design
reuses directly.

## Approach Options

Recommended: add `src/algorithm/redundancy.jl` with an internal
`SteinbergOptimizationCertificate`, a private construction helper, exact
factor-product replay, rewrite-log count validation, and a private verifier.
This follows the issue's suggested file boundary, keeps optimization concerns
out of `src/core/elementary_matrices.jl`, and creates the audit surface needed
before later shortening rules are added.

Alternative: place the certificate in `src/core/elementary_matrices.jl`. That
would keep it near canonical factor records, but the certificate is algorithmic
rewrite evidence rather than matrix construction behavior.

Alternative: export a public certificate constructor now. That conflicts with
the issue's internal-only scope and would prematurely freeze an optimizer API
before #188 decides public behavior.

## Chosen Design

Add a non-exported `SteinbergOptimizationCertificate` struct with these fields:

- `original_factors::Vector`
- `optimized_factors::Vector`
- `applied_rewrites::Vector`
- `comparison_summary`
- `original_product`
- `optimized_product`
- `verification`

Add `_steinberg_optimization_certificate(original_factors, optimized_factors,
applied_rewrite_metadata = ())` as the internal constructor. It collects the
input sequences, validates that the original sequence is nonempty, accepts only
ordinary-polynomial elementary or identity factors recognized by
`_canonical_elementary_factor_record`, requires every factor in both sequences
to share one size and base ring, and computes exact products by multiplying
matrix factors over that ring. Product inequality is not a construction error;
it is recorded as a failed verification status.

The comparison summary is a named tuple containing:

- `original_factor_count`, `optimized_factor_count`, and
  `factor_count_delta`;
- `original_metrics` and `optimized_metrics`, each with
  `max_elementary_factor_monomial_degree` and
  `total_elementary_factor_offdiagonal_monomials`;
- normalized `applied_rewrites`;
- `original_product`, `optimized_product`, `products_equal`, and
  `verification_status`.

Rewrite metadata stays generic but verifier-friendly. Each nonempty rewrite
record must expose a symbolic `rule_name`, an `original_factor_count`, and an
`optimized_factor_count`. The verifier checks that the sum of rewrite count
deltas equals the certificate's factor-count delta. This makes a no-op
certificate valid with an empty rule log and makes tampered logs with stale
counts fail verification.

Add `_verify_steinberg_optimization_certificate(certificate)::Bool`. It
recomputes sequence validation, exact products, metric summaries, the rewrite
count delta, and the expected comparison summary. It returns `false` for
replay/product/log/status mismatches and rethrows interrupts. Construction may
throw `ArgumentError` or `DimensionMismatch` for invalid matrices or mixed
rings before a certificate can be formed.

## Tests

Extend `test/expert/steinberg_factor_count_optimization.jl` in the existing
Steinberg testset:

- build a no-op certificate from a #288 catalog positive case and verify it;
- assert exact original and optimized products are equal;
- assert before/after factor counts are equal and `factor_count_delta == 0`;
- assert the summary records PR #179 metric values before and after
  optimization;
- assert the summary exposes the applied rewrite records, products, and true
  verification status;
- tamper one optimized factor over the same ring and require verification to
  return false;
- change the rule log count delta and require verification to return false;
- mix a factor from a different ordinary polynomial ring and require
  construction to throw `ArgumentError`.

Register no new public API tests because the certificate remains internal and
unexported.

## Verification

Focused verification:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Full verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must exit 0, including the no-op catalog certificate and
the three negative controls. The full package suite must remain green.

## Out Of Scope

Do not implement factor shortening. Do not export the certificate or helper
functions. Do not change `elementary_factorization(A)` or any public route to
optimize by default. Do not broaden Laurent or ToricBuilder support.

## Self-Review

This design has no placeholders, keeps construction internal, and ties every
issue requirement to either the certificate fields, comparison summary, exact
replay verifier, or focused expert tests. Product equality is checked by exact
matrix multiplication, while factor-count and rewrite-log checks remain
additional evidence rather than substitutes for equality.
