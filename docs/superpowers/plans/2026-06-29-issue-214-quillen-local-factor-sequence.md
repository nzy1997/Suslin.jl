# Issue 214 Quillen Local Factor Sequence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal replayable Quillen local factor sequence certificate that records ordered localized elementary factors with denominator and provenance metadata.

**Architecture:** Extend `src/algorithm/quillen_induction.jl` beside the existing single-correction local certificate. The new sequence certificate normalizes structured factor records into `QuillenLocalContribution` data, replays weighted elementary factors using the existing `_quillen_factors` path, and preserves denominator/provenance metadata for later patch assembly. Existing `QuillenLocalRealizationCertificate` values convert into length-one sequence certificates.

**Tech Stack:** Julia, Oscar, Suslin internal certificate structs, `Test`.

## Global Constraints

- Do not change the exported toy `QuillenPatch` API.
- Do not route through public `elementary_factorization`.
- Do not discover denominator covers, solve ideal membership, call the Murthy solver, or implement patch assembly.
- Do not implement the sequence certificate by merely storing already-materialized matrices.
- The verifier must replay elementary factors from `row`, `col`, `numerator`, and `denominator` data in order.
- The certificate must distinguish ordered factors, raw factor denominators, product denominator, replayed product/correction matrix, and normalized contribution/global factor data.
- Existing single-correction `QuillenLocalRealizationCertificate` behavior must keep working.
- Negative controls must reject corrupted factor entry, denominator, selected variable, factor provenance, and factor order.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add factor record, sequence certificate, replay verifier, constructor, and length-one conversion from existing local certificates.
- Create `test/expert/quillen_local_factor_sequence.jl`: two-factor #213 catalog coverage and negative controls.
- Modify `test/expert/quillen_local_certificate.jl`: add conversion regression for existing certificates.
- Modify `test/runtests.jl`: register the new expert file immediately before `expert/quillen_local_certificate.jl`.

### Task 1: Add Red Sequence-Certificate Tests

**Files:**
- Create: `test/expert/quillen_local_factor_sequence.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: future `Suslin.quillen_local_factor_sequence_certificate`, `Suslin.verify_quillen_local_factor_sequence_certificate`, `Suslin.QuillenLocalFactorSequenceCertificate`, and `Suslin.QuillenLocalElementaryFactor`.
- Produces: focused requirements for the implementation task.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/quillen_local_factor_sequence.jl` with helpers that load the #213 catalog, convert each catalog local factor into a structured record, rebuild tampered certificates, and check the positive and negative paths:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl"))

function qlfs_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function qlfs_factor_record(entry, index::Int)
    local_factor = entry.patch_case.local_factors[index]
    evidence = entry.local_evidence.records[index]
    return (;
        row = local_factor.correction.row,
        col = local_factor.correction.col,
        numerator = local_factor.correction.entry,
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        local_certificate = Suslin.LocalCertificate(
            local_factor.certificate.indices,
            local_factor.certificate.denominators,
        ),
        provenance = (;
            factor_index = index,
            local_index = index,
            fixture_id = entry.id,
            source = evidence.source,
            sequence_status = evidence.sequence_status,
        ),
        metadata = (;
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )
end
```

The positive test must build the sequence from
`QuillenMainlineFixtureCatalog.cases_by_id()["quillen-patched-substitution-witness-qq"]` and assert:

```julia
@test cert isa Suslin.QuillenLocalFactorSequenceCertificate
@test Suslin.verify_quillen_local_factor_sequence_certificate(cert)
@test length(cert.factors) == 2
@test cert.raw_denominators == [factor.denominator for factor in cert.factors]
@test cert.product_denominator == prod(cert.raw_denominators; init = one(cert.ring))
@test cert.normalized_global_elementary_factors == collect(entry.local_evidence.factors)
@test cert.local_product == qlfs_product(cert.normalized_global_elementary_factors, cert.ring, cert.size)
@test cert.local_correction == entry.local_evidence.expected_product
@test cert.verification.factor_provenance_ok
@test cert.verification.product_denominator_ok
@test cert.verification.overall_ok
```

The same file must rebuild the certificate with each of these corruptions and
assert verification returns `false`: changed first numerator, changed first
denominator, changed selected variable, changed first provenance record, and
reversed factor order.

- [ ] **Step 2: Register the red test**

In `test/runtests.jl`, add:

```julia
"expert/quillen_local_factor_sequence.jl",
```

immediately before:

```julia
"expert/quillen_local_certificate.jl",
```

- [ ] **Step 3: Run the focused red test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_local_factor_sequence.jl")'
```

Expected: FAIL with `UndefVarError` for `quillen_local_factor_sequence_certificate`.

### Task 2: Implement the Sequence Certificate

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: `_quillen_local_input_ring_size`, `_require_substitution_generator`, `_require_elementary_indices`, `_normalize_quillen_contribution`, `_quillen_factors`, `_quillen_product`, `_quillen_local_require_factor_matrix`, and `_quillen_local_patched_witness_summary`.
- Produces: `QuillenLocalElementaryFactor`, `QuillenLocalFactorSequenceCertificate`, `QuillenLocalFactorSequenceVerification`, `quillen_local_factor_sequence_certificate`, `replay_quillen_local_factor_sequence`, and `verify_quillen_local_factor_sequence_certificate`.

- [ ] **Step 1: Add the new structs**

Add the new structs after `QuillenLocalRealizationCertificate`:

```julia
struct QuillenLocalElementaryFactor
    row::Int
    col::Int
    numerator
    denominator
    coverage_multiplier
    provenance
    local_certificate::LocalCertificate
    metadata
end

struct QuillenLocalFactorSequenceVerification
    original_input
    selected_variable
    factor_count::Int
    raw_denominators::Vector
    product_denominator
    normalized_local_contributions::Vector{QuillenLocalContribution}
    normalized_global_elementary_factors::Vector
    local_product
    local_correction
    denominator_data::Vector{QuillenDenominatorData}
    factor_provenance::Vector
    factor_provenance_ok::Bool
    product_denominator_ok::Bool
    normalized_contributions_ok::Bool
    normalized_global_elementary_factors_ok::Bool
    local_product_ok::Bool
    local_correction_ok::Bool
    patched_substitution
    patched_substitution_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenLocalFactorSequenceCertificate
    original_input
    ring
    size::Int
    selected_variable
    factors::Vector{QuillenLocalElementaryFactor}
    raw_denominators::Vector
    product_denominator
    local_product
    local_correction
    normalized_local_contributions::Vector{QuillenLocalContribution}
    normalized_global_elementary_factors::Vector
    patched_substitution_witness
    chain_witness
    witness_metadata
    replay_metadata
    verification::QuillenLocalFactorSequenceVerification
end
```

- [ ] **Step 2: Add normalization helpers**

Add helper functions that require factor fields, coerce numerator and
denominator data into the certificate ring, derive `LocalCertificate([row, col],
[denominator, denominator])` when no local certificate is supplied, and validate
nonempty provenance:

```julia
function _quillen_local_sequence_factor_field(factor, field::Symbol)
    hasproperty(factor, field) ||
        throw(ArgumentError("local elementary factor missing field $(field)"))
    return getproperty(factor, field)
end
```

The helper `_quillen_local_sequence_factor(raw_factor, R, n::Int, index::Int)`
must return `QuillenLocalElementaryFactor(row, col, numerator, denominator,
coverage_multiplier, provenance, local_certificate, metadata)`.

- [ ] **Step 3: Add replay and constructor**

Add replay logic that:

```julia
contributions = [
    _normalize_quillen_contribution(
        QuillenLocalContribution(
            factor.local_certificate,
            factor.denominator,
            factor.coverage_multiplier,
            QuillenElementaryCorrection(factor.row, factor.col, factor.numerator),
        ),
        R,
        n,
    )
    for factor in certificate.factors
]
global_factors = _quillen_factors(R, n, contributions)
local_product = _quillen_product(R, n, global_factors)
```

The replay must recompute raw denominators, product denominator, denominator
data, ordered provenance, patched-substitution summary, replay metadata, and
all boolean checks. The verifier must return `false` on `ArgumentError` and
compare the stored verification with the replayed verification.

- [ ] **Step 4: Run the focused test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_local_factor_sequence.jl")'
```

Expected: PASS.

### Task 3: Add Length-One Conversion Coverage

**Files:**
- Modify: `test/expert/quillen_local_certificate.jl`

**Interfaces:**
- Consumes: `quillen_local_factor_sequence_certificate(certificate::QuillenLocalRealizationCertificate)`.
- Produces: regression coverage that existing local certificates still work and convert to one-factor sequences.

- [ ] **Step 1: Add the conversion assertion**

After the existing `correction_cert` assertions, add:

```julia
sequence_cert = Suslin.quillen_local_factor_sequence_certificate(correction_cert)
@test sequence_cert isa Suslin.QuillenLocalFactorSequenceCertificate
@test Suslin.verify_quillen_local_factor_sequence_certificate(sequence_cert)
@test length(sequence_cert.factors) == 1
@test sequence_cert.raw_denominators == [correction_cert.denominator]
@test sequence_cert.product_denominator == correction_cert.denominator
@test sequence_cert.local_product == correction_cert.local_product
@test sequence_cert.local_correction == correction_cert.local_correction
@test sequence_cert.normalized_global_elementary_factors == correction_cert.factors
@test sequence_cert.verification.overall_ok
```

- [ ] **Step 2: Run both required focused tests**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_local_factor_sequence.jl")'
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Expected: both commands exit 0.

### Task 4: Final Verification and Commit

**Files:**
- Modify: all implementation and test files changed by Tasks 1-3.

**Interfaces:**
- Consumes: completed implementation.
- Produces: committed worker-branch change ready for pull request creation.

- [ ] **Step 1: Run package entry-point verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 2: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: no output and exit 0.

- [ ] **Step 3: Commit implementation**

Run:

```bash
git add src/algorithm/quillen_induction.jl test/expert/quillen_local_factor_sequence.jl test/expert/quillen_local_certificate.jl test/runtests.jl docs/superpowers/plans/2026-06-29-issue-214-quillen-local-factor-sequence.md
git commit -m "feat: add quillen local factor sequence certificates"
```

## Final Verification

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_local_factor_sequence.jl")'
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: all commands exit 0.

## Self-Review

- Spec coverage: Tasks 1-3 cover ordered factors, raw denominators, product
  denominator, replayed product/correction, normalized contribution data,
  conversion, and negative controls.
- Placeholder scan: no `TBD`, `TODO`, unresolved fields, or incomplete steps are
  present.
- Type consistency: the plan uses the names
  `QuillenLocalFactorSequenceCertificate`,
  `QuillenLocalFactorSequenceVerification`,
  `QuillenLocalElementaryFactor`,
  `quillen_local_factor_sequence_certificate`,
  `replay_quillen_local_factor_sequence`, and
  `verify_quillen_local_factor_sequence_certificate` consistently.

## Execution Choice

Plan complete and saved to `docs/superpowers/plans/2026-06-29-issue-214-quillen-local-factor-sequence.md`. Two execution options:

1. Subagent-Driven (recommended) - dispatch a fresh subagent per task, review between tasks, fast iteration.
2. Inline Execution - execute tasks in this session using executing-plans, batch execution with checkpoints.

Automatic choice for this non-interactive run: Subagent-Driven, because it is marked recommended.
