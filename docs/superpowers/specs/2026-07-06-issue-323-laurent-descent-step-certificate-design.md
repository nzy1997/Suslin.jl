# Issue 323 Laurent Descent-Step Certificate Design

## Goal

Add an expert-only replayable Laurent descent-step certificate shell. A
certificate records the validated input profile and measure, one bounded-search
elementary Laurent operation, the replayed output profile and measure, and the
proof flag that the measure strictly decreases.

## Context

No repository `AGENTS.md` file is present in this worktree. Issue #323 has no
comments. Relevant merged PR context:

- #318 added the checked-in `case_008 d=14` Laurent boundary fixture.
- #320 added the d14 Laurent descent profile.
- #324 added the d14 Laurent descent measure contract.
- #325 added the d14 bounded Laurent elementary move search report.

Baseline verification in this sandbox:

- `julia --project=. -e 'using Pkg; Pkg.test()'` passed with public 729/729
  and internal 771/771.
- The #322 search report currently returns `status = :candidate_found` with
  225 replay-verified candidates. The first candidate is an entry addition
  from source 2 to target 1 with exponent `(-1, 1)`.

## Approach Options

Recommended: create `test/expert/laurent_descent_step_certificate.jl` as a
test-only certificate shell. Reuse the #322 elementary operation schema and
replay helper, recompute before and after profiles/measures from the actual
column, and expose a validator that returns `:ok` only when replay, stale-data
checks, and strict measure decrease all agree.

Alternative: move the certificate shell into `src/` as a production API. That
is premature for this issue because Laurent link witnesses, endpoint
reductions, and recursive peel integration are explicitly out of scope.

Alternative: only certify the real `case_008 d=14` search candidate. That would
cover the current fixture but would make the certificate shell harder to read
and slower to exercise. A small synthetic fixture should remain the primary
calibration case.

Chosen approach: expert-only test shell with a synthetic known-good case plus a
real `case_008 d=14` certificate constructor using the first #322 candidate.
This keeps the certificate shape stable while avoiding production API churn.

## Certificate Shape

Expose:

`laurent_descent_step_certificate(column, R, before_profile, operation) -> NamedTuple`

The returned named tuple has stable fields:

- `case_id`;
- `dimension`;
- `ring_generators`;
- `operation`;
- `before_profile`;
- `after_profile`;
- `before_measure`;
- `after_measure`;
- `status = :descent_step_certificate`;
- `replay_status = :ok`;
- `measure_relation = :strict_decrease`.

Also expose:

`case008_d14_laurent_descent_step_certificate() -> NamedTuple`

The `operation` field uses the same concrete schema as the #322 search report:

- `family = :entry_addition`;
- `target_index`;
- `source_index`;
- `coefficient`;
- `exponent = (a, b)`;
- `ring_generators = ("u", "v")`.

The certificate stores `ring_generators` both at the certificate level and on
the operation so replay can interpret the exponent tuple without guessing the
variable order.

## Replay Validation

Expose:

`validate_laurent_descent_step_certificate(cert, column, R)::Symbol`

The validator returns `:ok` only after:

- checking certificate-level fields, status, replay status, and measure
  relation;
- confirming the certificate and operation ring-generator metadata matches
  `gens(R)`;
- recomputing the before profile and before measure from `column`;
- replaying the recorded elementary operation against `column`;
- recomputing the after profile and after measure from the replayed column;
- comparing the recomputed after profile/measure with the certificate payload;
- proving `strictly_decreases_laurent_measure(before_measure, after_measure)`.

The validator never accepts an externally supplied after-profile or
after-measure without replaying the operation and recomputing both values.

## Profiles and Measures

For synthetic certificate tests, use a compact generic Laurent profile builder
with the same profile fields as the d14 profile shell:

- `case_id`;
- `dimension`;
- `ring_generators`;
- `nonzero_entries`;
- `max_entry_terms`;
- `entry_term_counts`;
- `valuation_ranges`;
- `newton_support_summary`;
- `leading_monomial_candidates`;
- `candidate_measure_families`;
- `status = :profile_only`.

Measures use the #321 component order and #322 support-set measure helper:

- `whole_support_count`;
- `max_entry_terms`;
- `valuation_span`;
- `leading_exponent`;
- `leading_entry_index`.

## Controls

The focused certificate test must cover:

- known-good synthetic descent certificate;
- tampered operation rejected after replay;
- stale after-profile rejected;
- equal-measure step rejected;
- wrong-ring certificate rejected;
- malformed operation schema rejected;
- certificate status changed away from `:descent_step_certificate` rejected;
- real `case_008 d=14` certificate from the first #322 candidate validates.

## Runner Registration

Register the new expert certificate file in `test/runtests.jl` immediately
after `expert/case008_d14_laurent_elementary_move_search.jl`.

## Verification

Required focused command:

```bash
julia --project=. -e 'include("test/expert/laurent_descent_step_certificate.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Also run:

```bash
git diff --check
```
