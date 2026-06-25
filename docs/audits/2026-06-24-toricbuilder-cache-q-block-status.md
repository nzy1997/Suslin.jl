# ToricBuilder Cache Q-Block Status Report

Date: 2026-06-25
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

| Case | Matrix size | Sparse nnz | Test level | Route status | Public elementary status | Determinant class | Decomposed base matrices | Runtime seconds |
| --- | ---: | ---: | --- | --- | --- | --- | ---: | ---: |
| case_001 | 6x6 | 30 | default_contract | gl_certificate_pass | staged_boundary | laurent_monomial_unit | 50 | 4.200 |
| case_002 | 14x14 | 79 | default_contract | gl_certificate_pass | staged_boundary | laurent_monomial_unit | 186 | 4.149 |
| case_003 | 6x6 | 27 | default_contract | gl_certificate_pass | staged_boundary | laurent_monomial_unit | 49 | 0.057 |
| case_004 | 18x18 | 73 | default_contract | gl_certificate_pass | staged_boundary | laurent_monomial_unit | 189 | 7.565 |
| case_005 | 14x14 | 90 | default_contract | gl_certificate_pass | staged_boundary | laurent_monomial_unit | 168 | 4.416 |
| case_006 | 18x18 | 99 | default_contract | gl_certificate_pass | staged_boundary | laurent_monomial_unit | 212 | 8.377 |
| case_007 | 42x42 | 546 | default_contract | not_exercised_in_default_report | not_run | not_run | not_run | not_run |
| case_008 | 30x30 | 477 | default_contract | not_exercised_in_default_report | not_run | not_run | not_run | not_run |
| case_009 | 62x62 | 739 | default_contract | not_exercised_in_default_report | not_run | not_run | not_run | not_run |
| case_010 | 6x6 | 34 | default_contract | gl_certificate_pass | staged_boundary | laurent_monomial_unit | 48 | 0.154 |
| case_011 | 288x288 | 14713 | optional_slow | not_exercised_in_default_report | not_run | not_run | not_run | not_run |
| case_012 | 98x98 | 1883 | default_contract | not_exercised_in_default_report | not_run | not_run | not_run | not_run |

## Exercised Evidence

- `case_001`: determinant `u*v`, route `gl_certificate_pass`, public `staged_boundary`, decomposed base matrices `50`, verified `true`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_002`: determinant `u*v^2`, route `gl_certificate_pass`, public `staged_boundary`, decomposed base matrices `186`, verified `true`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_003`: determinant `u*v`, route `gl_certificate_pass`, public `staged_boundary`, decomposed base matrices `49`, verified `true`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_004`: determinant `u^2*v^2`, route `gl_certificate_pass`, public `staged_boundary`, decomposed base matrices `189`, verified `true`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_005`: determinant `u^2*v`, route `gl_certificate_pass`, public `staged_boundary`, decomposed base matrices `168`, verified `true`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_006`: determinant `u*v^4`, route `gl_certificate_pass`, public `staged_boundary`, decomposed base matrices `212`, verified `true`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.
- `case_010`: determinant `u*v`, route `gl_certificate_pass`, public `staged_boundary`, decomposed base matrices `48`, verified `true`. normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is 1.

## Stage Timing Details

| Case | Determinant classification | Normalization | Certificate construction | Verification |
| --- | --- | --- | --- | --- |
| case_001 | pass (0.001s) | pass (0.002s) | pass (3.874s) | pass (0.014s) |
| case_002 | pass (0.007s) | pass (0.025s) | pass (3.139s) | pass (0.938s) |
| case_003 | pass (0.001s) | pass (0.003s) | pass (0.037s) | pass (0.013s) |
| case_004 | pass (0.009s) | pass (0.033s) | pass (5.459s) | pass (1.909s) |
| case_005 | pass (0.011s) | pass (0.159s) | pass (3.311s) | pass (0.893s) |
| case_006 | pass (0.023s) | pass (0.217s) | pass (5.771s) | pass (2.250s) |
| case_007 | not_run | not_run | not_run | not_run |
| case_008 | not_run | not_run | not_run | not_run |
| case_009 | not_run | not_run | not_run | not_run |
| case_010 | pass (0.001s) | pass (0.003s) | pass (0.129s) | pass (0.016s) |
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
