# Issue 327 Laurent Descent Internal Helpers Design

## Context

Issue #327 promotes reusable mechanics from the expert-only Laurent descent tests into internal algorithm code. The existing evidence lives in `test/expert/case008_d14_laurent_descent_measure_contract.jl`, `test/expert/case008_d14_laurent_elementary_move_search.jl`, and `test/expert/laurent_descent_step_certificate.jl`. The new code must stay case-agnostic and must not change the public API or claim full `case_008 d=14` support.

## Approach

Use a conservative internal-helper promotion. Add private helpers in `src/algorithm/column_reduction.jl` for two-generator Laurent columns:

- `_laurent_descent_measure_from_column(column, R; case_id = nothing)`
- `_strictly_decreases_laurent_measure(before, after)`
- `_replay_laurent_elementary_entry_addition(column, R, operation)`
- `_validate_laurent_descent_step_certificate(cert, column, R)`

These helpers compute the measure directly from a column's Laurent supports, replay the recorded elementary entry-addition operation over the supplied ring, and validate descent certificates by recomputing before and after measures from the input column and replayed operation. The implementation will accept optional diagnostic metadata such as `case_id`, but no reusable production helper will mention `case_008`, dimension 14, or fixture baseline constants.

## Alternatives Considered

1. Promote only the validator and leave measure/replay in tests. This keeps the production diff small but preserves duplicate mechanics and does not satisfy the issue's reusable-helper objective.
2. Promote the helper suite into `column_reduction.jl` without public exports. This satisfies the objective, preserves current API boundaries, and matches nearby internal ECP helper style. This is the chosen approach.
3. Create a new included source file for Laurent descent helpers. This would keep `column_reduction.jl` smaller, but the current module already keeps closely related ECP internals in that file and no include split is needed for this narrow change.

## Interfaces

The measure is a `NamedTuple` with stable fields:

- `status = :measure_contract`
- `order = :lexicographic_minimize`
- `components = (:whole_support_count, :max_entry_terms, :valuation_span, :leading_exponent, :leading_entry_index)`
- `whole_support_count`
- `max_entry_terms`
- `valuation_span`
- `leading_exponent`
- `leading_entry_index`
- `dimension`
- `ring_generators`
- `case_id` when supplied

The replay helper accepts an operation with fields `family`, `target_index`, `source_index`, `coefficient`, `exponent`, and `ring_generators`. It rejects missing generator metadata, malformed indices, equal target/source entries, wrong generator metadata, unsupported families, and exponent vectors whose length does not match the two-generator ring.

The validator accepts a certificate with stable fields matching the expert shell: `case_id`, `dimension`, `ring_generators`, `operation`, `before_measure`, `after_measure`, `status`, `replay_status`, and `measure_relation`. Optional `before_profile` and `after_profile` remain test-level metadata. Validation recomputes measures from the supplied column and the replayed after-column before comparing stored summaries.

## Test Design

Add `test/internal/laurent_descent_measure_helpers.jl`. The internal test will:

- compute the `case_008 d=14` fixture measure through the internal helper and assert the current baseline values;
- replay the known operation `(target_index = 1, source_index = 2, coefficient = 1, exponent = (-1, 1))`;
- assert that the after-measure strictly decreases;
- validate a minimal internal certificate built from replayed data;
- reject swapped ring generators, malformed source/target indices, target equal to source, stale supplied measures, zero-coefficient non-decreasing operation, and an operation whose after-measure is not recomputed from replay.

Update `test/runtests.jl` to register the new internal test. Update expert tests only where practical to call the internal helpers, preserving their existing negative controls and fixture-specific assertions.

## Out Of Scope

No public exports, no diagnostic support claim changes, and no implementation of Laurent link witnesses, endpoint reductions, determinant normalization, normality/conjugation replay, or recursive peel integration.

## Self Review

This design has no placeholders, keeps production helpers case-agnostic, preserves the public API, and maps every issue verification requirement to an internal test. The approach is narrow enough for one implementation plan.
