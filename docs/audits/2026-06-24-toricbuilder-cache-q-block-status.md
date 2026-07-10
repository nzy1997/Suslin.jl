# ToricBuilder Cache Q-Block Status Report

Date: 2026-06-29
Source fixture: `test/fixtures/toricbuilder_cache_q_blocks.jl`

This report records the current Suslin status for the checked-in ToricBuilder cache Q-block fixtures.
The default run fully exercises `case_001`, `case_002`, `case_003`, `case_004`, `case_005`, `case_006`, `case_010`; every other case is saved fixture data marked `not_exercised_in_default_report`.

## Summary

- Total checked-in Q-block cases: 12
- Fully exercised cases: 7
- GL certificate passes: 7
- Route errors: 0
- Not exercised in default report: 5

## Case Table

| Case | ToricBuilder a | ToricBuilder b | Matrix size | Sparse nnz | Test level | Status | Decomposed base matrices | Runtime seconds | Determinant | Max factor monomial degree | Total factor offdiag monomials |
| --- | --- | --- | ---: | ---: | --- | --- | ---: | ---: | --- | ---: | ---: |
| case_001 | x*y | x*y | 6x6 | 30 | default_contract | gl_certificate_pass / determinant_contract | 50 | 4.350 | u*v | 5 | 74 |
| case_002 | x^-1*y | x*y | 14x14 | 79 | default_contract | gl_certificate_pass / determinant_contract | 186 | 4.325 | u*v^2 | 7 | 281 |
| case_003 | x^2 | x^2 | 6x6 | 27 | default_contract | gl_certificate_pass / determinant_contract | 49 | 0.418 | u*v | 4 | 70 |
| case_004 | x^-1 | y^-1 | 18x18 | 73 | default_contract | gl_certificate_pass / determinant_contract | 189 | 7.660 | u^2*v^2 | 6 | 268 |
| case_005 | x*y | x*y^-1 | 14x14 | 90 | default_contract | gl_certificate_pass / determinant_contract | 168 | 4.498 | u^2*v | 5 | 257 |
| case_006 | x^-1 | x^3*y^2 | 18x18 | 99 | default_contract | gl_certificate_pass / determinant_contract | 212 | 8.166 | u*v^4 | 6 | 300 |
| case_007 | y^-2 | x^-2 | 42x42 | 546 | default_contract | not_exercised_in_default_report / not_run | not_run | not_run | not_run | not_run | not_run |
| case_008 | y^-2 | x^2 | 30x30 | 477 | default_contract | not_exercised_in_default_report / not_run | not_run | not_run | not_run | not_run | not_run |
| case_009 | x^-1*y | x^-1*y^-1 | 62x62 | 739 | default_contract | not_exercised_in_default_report / not_run | not_run | not_run | not_run | not_run | not_run |
| case_010 | x^-2*y^-1 | x^2*y | 6x6 | 34 | default_contract | gl_certificate_pass / determinant_contract | 48 | 0.134 | u*v | 5 | 95 |
| case_011 | x^-1*y^3 | x^3*y^-1 | 288x288 | 14713 | optional_slow | not_exercised_in_default_report / not_run | not_run | not_run | not_run | not_run | not_run |
| case_012 | x^-2 | x^-2*y^2 | 98x98 | 1883 | default_contract | not_exercised_in_default_report / not_run | not_run | not_run | not_run | not_run | not_run |

## Determinant Route Metadata

| Case | determinant_strategy | correction_side | determinant_source |
| --- | --- | --- | --- |
| case_001 | eager | not_run | not_run |
| case_002 | eager | not_run | not_run |
| case_003 | eager | not_run | not_run |
| case_004 | eager | not_run | not_run |
| case_005 | eager | not_run | not_run |
| case_006 | eager | not_run | not_run |
| case_007 | not_run | not_run | not_run |
| case_008 | not_run | not_run | not_run |
| case_009 | not_run | not_run | not_run |
| case_010 | eager | not_run | not_run |
| case_011 | not_run | not_run | not_run |
| case_012 | not_run | not_run | not_run |

## Exercised Evidence

- `case_001`: determinant `u*v`, route `gl_certificate_pass`, public `determinant_contract`, decomposed base matrices `50`, verified `true`, determinant_strategy `eager`, correction_side `not_run`, determinant_source `not_run`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_002`: determinant `u*v^2`, route `gl_certificate_pass`, public `determinant_contract`, decomposed base matrices `186`, verified `true`, determinant_strategy `eager`, correction_side `not_run`, determinant_source `not_run`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_003`: determinant `u*v`, route `gl_certificate_pass`, public `determinant_contract`, decomposed base matrices `49`, verified `true`, determinant_strategy `eager`, correction_side `not_run`, determinant_source `not_run`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_004`: determinant `u^2*v^2`, route `gl_certificate_pass`, public `determinant_contract`, decomposed base matrices `189`, verified `true`, determinant_strategy `eager`, correction_side `not_run`, determinant_source `not_run`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_005`: determinant `u^2*v`, route `gl_certificate_pass`, public `determinant_contract`, decomposed base matrices `168`, verified `true`, determinant_strategy `eager`, correction_side `not_run`, determinant_source `not_run`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_006`: determinant `u*v^4`, route `gl_certificate_pass`, public `determinant_contract`, decomposed base matrices `212`, verified `true`, determinant_strategy `eager`, correction_side `not_run`, determinant_source `not_run`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_010`: determinant `u*v`, route `gl_certificate_pass`, public `determinant_contract`, decomposed base matrices `48`, verified `true`, determinant_strategy `eager`, correction_side `not_run`, determinant_source `not_run`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.

## Stage Timing Details

| Case | Determinant classification | Normalization | Certificate construction | Verification |
| --- | --- | --- | --- | --- |
| case_001 | pass (0.001s) | pass (0.003s) | pass (3.947s) | pass (0.018s) |
| case_002 | pass (0.008s) | pass (0.125s) | pass (3.106s) | pass (1.045s) |
| case_003 | pass (0.001s) | pass (0.003s) | pass (0.396s) | pass (0.014s) |
| case_004 | pass (0.009s) | pass (0.036s) | pass (5.913s) | pass (1.667s) |
| case_005 | pass (0.012s) | pass (0.173s) | pass (3.746s) | pass (0.522s) |
| case_006 | pass (0.023s) | pass (0.220s) | pass (5.862s) | pass (1.949s) |
| case_007 | not_run | not_run | not_run | not_run |
| case_008 | not_run | not_run | not_run | not_run |
| case_009 | not_run | not_run | not_run | not_run |
| case_010 | pass (0.001s) | pass (0.003s) | pass (0.110s) | pass (0.016s) |
| case_011 | not_run | not_run | not_run | not_run |
| case_012 | not_run | not_run | not_run | not_run |

## Not Exercised Boundary

The non-exercised rows are real stored matrices, not placeholders. They are intentionally not routed through the full factorization stack in the default report because large Laurent GL certification can be much slower than fixture/schema validation.
Use `--exercise=case_001,case_003,...` to expand the set of cases that the report probes.

## Reproduction

```text
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl
julia --project=. test/internal/toricbuilder_cache_status_report.jl
```
