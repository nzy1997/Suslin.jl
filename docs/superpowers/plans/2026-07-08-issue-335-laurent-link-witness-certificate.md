# Issue 335 Laurent Link-Witness Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an expert-only replayable Laurent link-witness certificate shell.

**Architecture:** Keep the shell in `test/expert/laurent_link_witness_certificate.jl`, consuming the #334 search helper and reusing its witness replay helpers. The certificate constructor records one replayed witness, while the validator recomputes endpoints from the supplied source column and rejects stale or malformed certificate data.

**Tech Stack:** Julia, Test, Oscar, Suslin.jl expert test helpers.

## Global Constraints

- Do not add a production Laurent ECP stage or public API.
- Do not implement endpoint reductions, normality/conjugation replay, determinant normalization, recursive peel integration, or full `case_008` success.
- The certificate status is exactly `:link_witness_certificate`.
- The accepted certificate has `replay_status == :ok`, `identity_status == :verified`, and `next_boundary == :laurent_endpoint_reduction`.
- The validator signature is `validate_laurent_link_witness_certificate(cert, column, R)::Symbol`.
- The validator must recompute endpoint data from witness replay and must not trust supplied endpoint fields.
- The d14 summary status is `:link_witness_certificate` only when the bounded d14 search report found a replay-verified candidate that validates through this shell.
- The d14 summary status is `:no_link_witness_candidate` when the bounded d14 search is exhausted.
- Synthetic-only success must not claim real d14 endpoint-reduction progress.

---

### Task 1: Add the Failing Certificate Contract

**Files:**
- Create: `test/expert/laurent_link_witness_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `case008_d14_laurent_link_witness_search_report()`
- Consumes: `_laurent_link_witness_candidate_from_replay(column, R, witness; case_id, require_strict = true)`
- Consumes: `verify_laurent_link_witness_candidate(column, R, candidate)::Bool`
- Produces later: `laurent_link_witness_certificate(context, witness, column, R)::NamedTuple`
- Produces later: `validate_laurent_link_witness_certificate(cert, column, R)::Symbol`
- Produces later: `case008_d14_laurent_link_witness_certificate_summary()::NamedTuple`

- [ ] **Step 1: Create the failing test file**

Create `test/expert/laurent_link_witness_certificate.jl` with this initial content:

```julia
using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_link_witness_search.jl"))
end

@testset "Laurent link-witness certificate shell" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/laurent_link_witness_certificate.jl\"",
        runtests,
    )

    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [one(R), u + v]
    context = (;
        case_id = "synthetic",
        dimension = 2,
        ring_generators = ("u", "v"),
        status = :link_witness_context,
    )
    witness = (;
        family = :two_entry_laurent_combination,
        pivot_index = 2,
        partner_index = 1,
        coefficient = 1,
        exponent = (1, 0),
        ring_generators = ("u", "v"),
    )

    cert = laurent_link_witness_certificate(context, witness, column, R)
    @test cert.case_id == "synthetic"
    @test cert.dimension == 2
    @test cert.ring_generators == ("u", "v")
    @test cert.context_status == :link_witness_context
    @test cert.witness == witness
    @test cert.source_endpoint.leading_exponent == (1, 0)
    @test cert.target_endpoint.leading_exponent == (0, 1)
    @test cert.replay_status == :ok
    @test cert.identity_status == :verified
    @test cert.next_boundary == :laurent_endpoint_reduction
    @test cert.status == :link_witness_certificate
    @test validate_laurent_link_witness_certificate(cert, column, R) == :ok
end

@testset "case_008 d=14 Laurent link-witness certificate summary" begin
    summary = case008_d14_laurent_link_witness_certificate_summary()
    @test summary.case_id == "case_008"
    @test summary.dimension == 14
    @test summary.search_status in (:candidate_found, :exhausted)
    @test summary.d14_status in (:link_witness_certificate, :no_link_witness_candidate)

    if summary.search_status == :candidate_found
        @test summary.d14_status == :link_witness_certificate
        @test summary.certificate.case_id == "case_008"
        @test validate_laurent_link_witness_certificate(
            summary.certificate,
            summary.source_column,
            summary.ring,
        ) == :ok
    else
        @test summary.d14_status == :no_link_witness_candidate
        @test !hasproperty(summary, :certificate)
    end
end
```

- [ ] **Step 2: Register the file**

In `test/runtests.jl`, insert this line immediately after
`"expert/case008_d14_laurent_link_witness_search.jl",`:

```julia
        "expert/laurent_link_witness_certificate.jl",
```

- [ ] **Step 3: Run the focused test and confirm the expected failure**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_link_witness_certificate.jl")'
```

Expected: nonzero exit with `UndefVarError: laurent_link_witness_certificate not defined`.

---

### Task 2: Implement Certificate Construction and Validation

**Files:**
- Modify: `test/expert/laurent_link_witness_certificate.jl`

**Interfaces:**
- Consumes: Task 1 test contract.
- Produces: `laurent_link_witness_certificate(context, witness, column, R)::NamedTuple`
- Produces: `validate_laurent_link_witness_certificate(cert, column, R)::Symbol`
- Produces: `_laurent_link_witness_certificate_from_replay(context, witness, column, R; require_strict = true)::NamedTuple`

- [ ] **Step 1: Add certificate constants and field helpers**

Add these definitions after the include block:

```julia
const LAURENT_LINK_WITNESS_CERTIFICATE_STATUS = :link_witness_certificate
const LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY =
    :laurent_endpoint_reduction

const LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :context_status,
    :witness,
    :source_endpoint,
    :target_endpoint,
    :replay_status,
    :identity_status,
    :next_boundary,
    :status,
)

function _laurent_link_witness_certificate_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end
```

- [ ] **Step 2: Implement the replay-derived certificate constructor**

Add this implementation below the constants:

```julia
function _laurent_link_witness_certificate_from_replay(
    context,
    witness,
    column,
    R;
    require_strict::Bool = true,
)
    _laurent_link_witness_certificate_has_fields(
        context,
        (:case_id, :dimension, :ring_generators, :status),
    ) || throw(ArgumentError("context must include case_id, dimension, ring_generators, and status"))
    context.status == :link_witness_context ||
        throw(ArgumentError("context must have status :link_witness_context"))
    context.dimension == length(column) ||
        throw(ArgumentError("context dimension does not match source column"))
    context.ring_generators == _ring_generator_names(R) ||
        throw(ArgumentError("context ring_generators do not match the ring"))

    candidate = _laurent_link_witness_candidate_from_replay(
        column,
        R,
        witness;
        case_id = context.case_id,
        require_strict,
    )
    verify_laurent_link_witness_candidate(column, R, candidate) ||
        throw(ArgumentError("witness replay did not verify as a link witness"))

    return (;
        case_id = context.case_id,
        dimension = context.dimension,
        ring_generators = context.ring_generators,
        context_status = context.status,
        witness = candidate.witness,
        source_endpoint = candidate.source_endpoint,
        target_endpoint = candidate.target_endpoint,
        replay_status = :ok,
        identity_status = :verified,
        next_boundary = LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY,
        status = LAURENT_LINK_WITNESS_CERTIFICATE_STATUS,
    )
end

function laurent_link_witness_certificate(context, witness, column, R)
    return _laurent_link_witness_certificate_from_replay(
        context,
        witness,
        column,
        R;
        require_strict = true,
    )
end
```

- [ ] **Step 3: Implement the validator**

Add this implementation below the constructor:

```julia
function validate_laurent_link_witness_certificate(cert, column, R)::Symbol
    try
        _laurent_link_witness_certificate_has_fields(
            cert,
            LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS,
        ) || return :missing_certificate_fields
        cert.status == LAURENT_LINK_WITNESS_CERTIFICATE_STATUS ||
            return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.identity_status == :verified || return :wrong_identity_status
        cert.next_boundary == LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY ||
            return :wrong_next_boundary
        cert.dimension == length(column) || return :wrong_dimension
        cert.ring_generators == _ring_generator_names(R) ||
            return :wrong_ring_generators
        cert.context_status == :link_witness_context ||
            return :wrong_context_status

        context = (;
            case_id = cert.case_id,
            dimension = cert.dimension,
            ring_generators = cert.ring_generators,
            status = cert.context_status,
        )
        candidate = _laurent_link_witness_candidate_from_replay(
            column,
            R,
            cert.witness;
            case_id = cert.case_id,
            require_strict = true,
        )
        candidate.source_endpoint == cert.source_endpoint ||
            return :stale_source_endpoint
        candidate.target_endpoint == cert.target_endpoint ||
            return :stale_target_endpoint
        verify_laurent_link_witness_candidate(column, R, candidate) ||
            return :identity_replay_failed

        expected = _laurent_link_witness_certificate_from_replay(
            context,
            cert.witness,
            column,
            R;
            require_strict = true,
        )
        cert == expected || return :stale_certificate
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :identity_replay_failed
    end
end
```

- [ ] **Step 4: Run the focused test and confirm it passes**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_link_witness_certificate.jl")'
```

Expected: the first certificate testset passes; the d14 summary test may still
fail with `UndefVarError: case008_d14_laurent_link_witness_certificate_summary not defined`.

---

### Task 3: Add Negative Controls and d14 Consumption

**Files:**
- Modify: `test/expert/laurent_link_witness_certificate.jl`

**Interfaces:**
- Consumes: Task 2 certificate constructor and validator.
- Produces: `case008_d14_laurent_link_witness_certificate_summary(fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture())::NamedTuple`

- [ ] **Step 1: Add the d14 summary helper**

Add this function below the certificate validator:

```julia
function case008_d14_laurent_link_witness_certificate_summary(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    source = _case008_d14_link_witness_source_data(fixture)
    search_report = case008_d14_laurent_link_witness_search_report(fixture)
    search_validation =
        validate_case008_d14_laurent_link_witness_search_report(
            search_report,
            fixture,
        )
    search_validation == :ok ||
        throw(ArgumentError("d14 link-witness search report must validate; got $(search_validation)"))

    common = (;
        case_id = search_report.case_id,
        dimension = search_report.dimension,
        search_status = search_report.status,
        search_next_boundary = search_report.next_boundary,
        candidate_count = search_report.candidate_count,
        replay_verified_count = search_report.replay_verified_count,
        source_column = source.replay.after_column,
        ring = fixture.ring,
    )

    if search_report.status == :candidate_found
        candidate = first(search_report.candidates)
        cert = laurent_link_witness_certificate(
            source.context,
            candidate.witness,
            source.replay.after_column,
            fixture.ring,
        )
        validate_laurent_link_witness_certificate(
            cert,
            source.replay.after_column,
            fixture.ring,
        ) == :ok ||
            throw(ArgumentError("first d14 candidate did not validate through the certificate shell"))
        return merge(
            common,
            (;
                d14_status = :link_witness_certificate,
                certificate = cert,
            ),
        )
    end

    return merge(common, (; d14_status = :no_link_witness_candidate))
end
```

- [ ] **Step 2: Add negative-control assertions**

Inside the `Laurent link-witness certificate shell` testset, after the positive
validator assertion, add:

```julia
    tampered_coefficient = merge(
        cert,
        (; witness = merge(witness, (; coefficient = 0))),
    )
    @test validate_laurent_link_witness_certificate(
        tampered_coefficient,
        column,
        R,
    ) == :identity_replay_failed

    wrong_generators = merge(cert, (; ring_generators = ("v", "u")))
    @test validate_laurent_link_witness_certificate(
        wrong_generators,
        column,
        R,
    ) == :wrong_ring_generators

    stale_source = merge(
        cert,
        (;
            source_endpoint = merge(
                cert.source_endpoint,
                (; term_count = cert.source_endpoint.term_count + 1),
            ),
        ),
    )
    @test validate_laurent_link_witness_certificate(
        stale_source,
        column,
        R,
    ) == :stale_source_endpoint

    stale_target = merge(
        cert,
        (;
            target_endpoint = merge(
                cert.target_endpoint,
                (; term_count = cert.target_endpoint.term_count + 1),
            ),
        ),
    )
    @test validate_laurent_link_witness_certificate(
        stale_target,
        column,
        R,
    ) == :stale_target_endpoint

    malformed_witness = merge(
        cert,
        (; witness = _laurent_link_without_field(witness, :family)),
    )
    @test validate_laurent_link_witness_certificate(
        malformed_witness,
        column,
        R,
    ) == :identity_replay_failed

    nonidentity_cert = _laurent_link_witness_certificate_from_replay(
        context,
        merge(witness, (; exponent = (0, 0))),
        column,
        R;
        require_strict = false,
    )
    @test validate_laurent_link_witness_certificate(
        nonidentity_cert,
        column,
        R,
    ) == :identity_replay_failed
```

- [ ] **Step 3: Tighten the d14 summary assertions**

In the d14 testset, add:

```julia
    if summary.search_status == :candidate_found
        @test summary.search_next_boundary == :laurent_link_witness_certificate
        @test summary.candidate_count > 0
        @test summary.replay_verified_count == summary.candidate_count
        @test summary.certificate.next_boundary == :laurent_endpoint_reduction
        @test summary.certificate.status == :link_witness_certificate
    else
        @test summary.search_next_boundary == :laurent_link_witness_search_expansion
        @test summary.candidate_count == 0
        @test summary.replay_verified_count == 0
    end
```

- [ ] **Step 4: Run the focused test and confirm all new behavior passes**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_link_witness_certificate.jl")'
```

Expected: command exits 0.

- [ ] **Step 5: Run the package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

---

### Task 4: Final Verification and Commit

**Files:**
- Modify: `test/expert/laurent_link_witness_certificate.jl`
- Modify: `test/runtests.jl`
- Modify: `docs/superpowers/plans/2026-07-08-issue-335-laurent-link-witness-certificate.md`

**Interfaces:**
- Consumes: Tasks 1-3.
- Produces: committed implementation ready for PR.

- [ ] **Step 1: Run formatting and whitespace checks**

Run:

```bash
git diff --check
```

Expected: command exits 0.

- [ ] **Step 2: Review the final diff**

Run:

```bash
git diff --stat
git diff -- test/expert/laurent_link_witness_certificate.jl test/runtests.jl docs/superpowers/plans/2026-07-08-issue-335-laurent-link-witness-certificate.md
```

Expected: the diff only contains the certificate expert test, test registration,
and this plan.

- [ ] **Step 3: Commit the implementation**

Run:

```bash
git add test/expert/laurent_link_witness_certificate.jl test/runtests.jl docs/superpowers/plans/2026-07-08-issue-335-laurent-link-witness-certificate.md
git commit -m "Add issue 335 Laurent link witness certificate shell"
```

Expected: commit succeeds.

## Self-Review

- Spec coverage: Tasks 1-3 cover the stable certificate fields, replay-derived
  validator, synthetic positive certificate, negative controls, d14 conditional
  summary, and test registration.
- Placeholder scan: no deferred implementation placeholders are present.
- Type consistency: all produced interfaces use NamedTuple values and the
  validator returns `Symbol`.
