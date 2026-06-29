# Issue 215 Quillen Denominator Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a checked internal extraction stage that turns verified Quillen local factor-sequence certificates into raw denominator-cover candidates.

**Architecture:** Extend `src/algorithm/quillen_induction.jl` beside the #214 factor-sequence verifier. The new candidate record stores per-local product support denominators, exact factor entries and provenance, and replay metadata; verification recomputes the candidate from the stored local sequence certificates.

**Tech Stack:** Julia, Oscar, Suslin internal Quillen certificate structs, `Test`.

## Global Constraints

- Do not prove that extracted denominators cover the ring.
- Do not choose an exponent `l` or verify `sum g_i * r_i^l == 1`.
- Do not assemble global factors or call Murthy local solving.
- Read structured `QuillenLocalFactorSequenceCertificate` records, not materialized matrices, when extracting denominators.
- Keep raw denominators `r_i` separate from later cover powers.
- Reject empty input, unverified local sequence certificates, mixed original input, mixed ring, mixed size, and mixed selected variable.
- Record exact factor-level denominator and provenance data so verifier replay rejects manual edits after extraction.
- For #215, local support denominators use the #214 product denominator and record support kind `:product`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add support/candidate structs, extraction, replay, equality helpers, and verification.
- Create `test/expert/quillen_denominator_extraction.jl`: positive two-open extraction from #213/#214 evidence and negative controls.
- Modify `test/runtests.jl`: register the new expert file after `expert/quillen_local_factor_sequence.jl`.

### Task 1: Add Red Denominator-Extraction Tests

**Files:**
- Create: `test/expert/quillen_denominator_extraction.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: future `Suslin.extract_quillen_denominator_cover_candidate`, `Suslin.verify_quillen_denominator_cover_candidate`, `Suslin.replay_quillen_denominator_cover_candidate`, `Suslin.QuillenDenominatorCoverCandidate`, and `Suslin.QuillenLocalDenominatorSupport`.
- Produces: failing behavioral coverage for the implementation task.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/quillen_denominator_extraction.jl` with helpers that load
`test/fixtures/quillen_mainline_cases.jl`, build one verified #214 sequence
certificate per local factor from the `quillen-two-open-cover-qq` entry, and
exercise the extractor.

The helper for factor records must mirror the #214 record shape:

```julia
function qde_factor_record(entry, index::Int)
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
            factor_index = 1,
            sequence_index = 1,
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

The positive test must assert:

```julia
@test candidate isa Suslin.QuillenDenominatorCoverCandidate
@test Suslin.verify_quillen_denominator_cover_candidate(candidate)
@test candidate.raw_denominators == [data.denominator for data in entry.patch_case.denominator_data]
@test candidate.selected_variable == entry.patch_case.substitution_variable
@test candidate.original_input == entry.patch_case.target_matrix
@test length(candidate.local_supports) == 2
@test all(support -> support.support_kind == :product, candidate.local_supports)
@test all(support -> support.replay_ok, candidate.local_supports)
@test [only(support.factor_denominators) for support in candidate.local_supports] == candidate.raw_denominators
@test candidate.verification.local_certificates_ok
@test candidate.verification.local_supports_ok
@test candidate.verification.raw_denominators_ok
@test candidate.verification.overall_ok
```

The negative controls must assert construction or verification rejects:

```julia
@test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(Suslin.QuillenLocalFactorSequenceCertificate[])
@test !Suslin.verify_quillen_denominator_cover_candidate(dropped_candidate)
@test !Suslin.verify_quillen_denominator_cover_candidate(edited_raw_denominator_candidate)
@test !Suslin.verify_quillen_denominator_cover_candidate(edited_support_denominator_candidate)
@test !Suslin.verify_quillen_denominator_cover_candidate(edited_factor_denominator_candidate)
@test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(mixed_variable_certificates)
```

- [ ] **Step 2: Register the red test**

In `test/runtests.jl`, add:

```julia
"expert/quillen_denominator_extraction.jl",
```

immediately after:

```julia
"expert/quillen_local_factor_sequence.jl",
```

- [ ] **Step 3: Run the focused red test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_extraction.jl")'
```

Expected: FAIL with `UndefVarError` for `QuillenDenominatorCoverCandidate` or
`extract_quillen_denominator_cover_candidate`.

### Task 2: Implement Denominator-Cover Candidate Extraction

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: `QuillenLocalFactorSequenceCertificate`, `QuillenLocalElementaryFactor`, `verify_quillen_local_factor_sequence_certificate`, `_require_supported_quillen_ring`, `_require_substitution_generator`, `_quillen_local_sequence_original_input`, `_quillen_local_sequence_factor`, and `_same_quillen_local_factor_sequence_verification`.
- Produces: `QuillenLocalDenominatorSupport`, `QuillenDenominatorCoverCandidateVerification`, `QuillenDenominatorCoverCandidate`, `extract_quillen_denominator_cover_candidate`, `replay_quillen_denominator_cover_candidate`, and `verify_quillen_denominator_cover_candidate`.

- [ ] **Step 1: Add candidate structs**

Add the structs after `QuillenLocalFactorSequenceCertificate`:

```julia
struct QuillenLocalDenominatorSupport
    local_index::Int
    support_denominator
    support_kind::Symbol
    factor_denominators::Vector
    factor_entries::Vector{QuillenLocalElementaryFactor}
    factor_provenance::Vector
    replayed_denominator
    replay_equality
    replay_ok::Bool
end

struct QuillenDenominatorCoverCandidateVerification
    local_count::Int
    raw_denominators::Vector
    local_certificates_ok::Bool
    same_original_input_ok::Bool
    same_ring_ok::Bool
    same_size_ok::Bool
    same_selected_variable_ok::Bool
    local_supports::Vector{QuillenLocalDenominatorSupport}
    local_supports_ok::Bool
    raw_denominators_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenDenominatorCoverCandidate
    original_input
    ring
    size::Int
    selected_variable
    local_certificates::Vector{QuillenLocalFactorSequenceCertificate}
    raw_denominators::Vector
    local_supports::Vector{QuillenLocalDenominatorSupport}
    replay_metadata
    verification::QuillenDenominatorCoverCandidateVerification
end
```

- [ ] **Step 2: Add support replay helpers**

Implement `_quillen_local_denominator_support(certificate, local_index::Int)`
so it normalizes the certificate factors with `_quillen_local_sequence_factor`,
copies exact factor entries and provenance, computes:

```julia
replayed_denominator = prod(factor.denominator for factor in factors; init = one(certificate.ring))
replay_ok = certificate.raw_denominators == factor_denominators &&
            certificate.product_denominator == replayed_denominator
```

and returns a `QuillenLocalDenominatorSupport` with `support_kind = :product`
and `support_denominator = certificate.product_denominator`.

Add `_same_quillen_local_denominator_support` and
`_same_quillen_local_denominator_supports` helpers that compare all stored
fields, including `factor_entries`.

- [ ] **Step 3: Add candidate replay and verifier**

Implement `replay_quillen_denominator_cover_candidate(candidate)` to recompute
context alignment, support records, raw denominators, and replay metadata from
`candidate.local_certificates`. It must return a
`QuillenDenominatorCoverCandidateVerification` whose `overall_ok` requires all
local certificates to verify, all context checks to pass, every local support
to replay, stored supports to match replayed supports, stored raw denominators
to match replayed raw denominators, and replay metadata to match.

Implement `verify_quillen_denominator_cover_candidate(candidate)::Bool` so it
returns `false` on `ArgumentError` and requires both `replay.overall_ok` and a
stored-vs-replayed verification match.

- [ ] **Step 4: Add extraction constructor**

Implement `extract_quillen_denominator_cover_candidate(certificates)` to reject
empty input, reject unverified certificates, reject mixed context, compute
supports and raw denominators, build replay metadata, replay the candidate, and
throw `ArgumentError("Quillen denominator-cover candidate data does not replay")`
if the replay does not pass.

- [ ] **Step 5: Run the focused test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_extraction.jl")'
```

Expected: PASS.

### Task 3: Final Verification, Review, and Commit

**Files:**
- Modify: all files changed by Tasks 1-2.

**Interfaces:**
- Consumes: completed extraction implementation.
- Produces: committed worker-branch change ready for pull request creation.

- [ ] **Step 1: Run focused issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_extraction.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package entry-point verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: no output and exit 0.

- [ ] **Step 4: Commit implementation**

Run:

```bash
git add src/algorithm/quillen_induction.jl test/expert/quillen_denominator_extraction.jl test/runtests.jl docs/superpowers/specs/2026-06-30-issue-215-quillen-denominator-extraction-design.md docs/superpowers/plans/2026-06-30-issue-215-quillen-denominator-extraction.md
git commit -m "feat: extract quillen denominator cover candidates"
```

## Final Verification

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_extraction.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: all commands exit 0.

## Self-Review

- Spec coverage: Tasks 1-2 cover nonempty verified local sequence input, context
  alignment, raw per-local denominators, support-kind recording, factor-level
  provenance, replay metadata, and negative controls.
- Placeholder scan: no placeholders or incomplete steps remain.
- Type consistency: the plan uses `QuillenLocalDenominatorSupport`,
  `QuillenDenominatorCoverCandidateVerification`,
  `QuillenDenominatorCoverCandidate`,
  `extract_quillen_denominator_cover_candidate`,
  `replay_quillen_denominator_cover_candidate`, and
  `verify_quillen_denominator_cover_candidate` consistently.

## Execution Choice

Plan complete and saved to `docs/superpowers/plans/2026-06-30-issue-215-quillen-denominator-extraction.md`. Two execution options:

1. Subagent-Driven (recommended) - dispatch a fresh subagent per task, review between tasks, fast iteration.
2. Inline Execution - execute tasks in this session using executing-plans, batch execution with checkpoints.

Automatic choice for this non-interactive run: Subagent-Driven, because it is marked recommended.
