# Issue 339 Laurent Link-Witness Internal Helpers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote reusable Laurent link-witness replay, candidate, endpoint, and certificate validation mechanics into internal Suslin helpers.

**Architecture:** Add generic internal helpers to `src/algorithm/column_reduction.jl` beside the existing `_laurent_descent_*` helpers. Keep d14 fixture constants and reports in expert tests, and update expert code to delegate generic witness/certificate mechanics to the new `Suslin._...` helpers.

**Tech Stack:** Julia, Oscar, Suslin.jl internal helpers, `Test`.

## Global Constraints

- Do not expose a public API.
- Do not claim full `case_008 d=14` support.
- Do not bake d14 fixture fingerprints, fixture constants, or `case008_d14_*` helper names into generic internal helper code.
- Do not update `diagnose_unimodular_column_reduction` except as needed to compile shared helpers.
- Do not implement Laurent endpoint reductions, normality/conjugation replay, determinant normalization, recursive peel integration, or full `case_008` support.
- Do not use third-party attachment patches or ZIP files as implementation sources.
- The supported witness family is exactly `:two_entry_laurent_combination`.
- Replay semantics are `pivot_entry + coefficient * u^a * v^b * partner_entry`, with `(a, b)` in `ring_generators` order.
- `_validate_laurent_link_witness_certificate` must recompute endpoint metadata from replay rather than trusting supplied endpoint fields.
- The internal test must assert the first d14 candidate has `case_id == "case_008"`, `pivot_index == 10`, `partner_index == 1`, `coefficient == 1`, `exponent == (1, -1)`, and `next_boundary == :laurent_endpoint_reduction`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add internal link-witness constants, status checks, replay/candidate helpers, endpoint metadata, and certificate validator.
- Create `test/internal/laurent_link_witness_helpers.jl`: focused synthetic and narrow d14 contract tests for the promoted helpers.
- Modify `test/runtests.jl`: register the new internal test immediately after `internal/laurent_descent_measure_helpers.jl`.
- Modify `test/expert/case008_d14_laurent_link_witness_search.jl`: keep d14 search/report code but forward generic helper names to `Suslin._...`.
- Modify `test/expert/laurent_link_witness_certificate.jl`: keep expert summary code but forward construction and validation to `Suslin._...`.

---

### Task 1: Add the Failing Internal Contract

**Files:**
- Create: `test/internal/laurent_link_witness_helpers.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin._laurent_link_witness_status(witness, n, R)::Symbol`
- Consumes: `Suslin._laurent_link_witness_candidate_from_replay(column, R, witness; case_id, require_strict = true)::NamedTuple`
- Consumes: `Suslin._verify_laurent_link_witness_candidate(column, R, candidate)::Bool`
- Consumes: `Suslin._laurent_link_witness_certificate_from_replay(context, witness, column, R; require_strict = true)::NamedTuple`
- Consumes: `Suslin._validate_laurent_link_witness_certificate(cert, column, R)::Symbol`

- [ ] **Step 1: Create the failing internal test**

Create `test/internal/laurent_link_witness_helpers.jl`:

```julia
using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_search_report)
    include(joinpath(@__DIR__, "..", "expert", "case008_d14_laurent_link_witness_search.jl"))
end

function _internal_link_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

@testset "internal Laurent link-witness helpers" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin("\"internal/laurent_link_witness_helpers.jl\"", runtests)

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

    @test Suslin._laurent_link_witness_status(witness, length(column), R) == :ok
    candidate = Suslin._laurent_link_witness_candidate_from_replay(
        column,
        R,
        witness;
        case_id = context.case_id,
    )
    @test Suslin._verify_laurent_link_witness_candidate(column, R, candidate)
    @test candidate.source_endpoint.leading_exponent == (1, 0)
    @test candidate.target_endpoint.leading_exponent == (0, 1)
    @test candidate.measure_relation == :strict_decrease

    cert = Suslin._laurent_link_witness_certificate_from_replay(
        context,
        witness,
        column,
        R,
    )
    @test cert.case_id == "synthetic"
    @test cert.next_boundary == :laurent_endpoint_reduction
    @test cert.status == :link_witness_certificate
    @test Suslin._validate_laurent_link_witness_certificate(cert, column, R) == :ok

    missing_field = _internal_link_without_field(witness, :family)
    @test Suslin._laurent_link_witness_status(missing_field, length(column), R) ==
          :malformed_witness
    missing_field_cert = merge(cert, (; witness = missing_field))
    @test Suslin._validate_laurent_link_witness_certificate(
        missing_field_cert,
        column,
        R,
    ) == :identity_replay_failed

    equal_indices = merge(witness, (; partner_index = 2))
    @test Suslin._laurent_link_witness_status(equal_indices, length(column), R) ==
          :pivot_partner_index_equality
    @test !Suslin._verify_laurent_link_witness_candidate(
        column,
        R,
        merge(candidate, (; witness = equal_indices)),
    )

    wrong_generators = merge(witness, (; ring_generators = ("v", "u")))
    @test Suslin._laurent_link_witness_status(wrong_generators, length(column), R) ==
          :wrong_ring_generators

    stale_source = merge(
        cert,
        (;
            source_endpoint = merge(
                cert.source_endpoint,
                (; term_count = cert.source_endpoint.term_count + 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_link_witness_certificate(
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
    @test Suslin._validate_laurent_link_witness_certificate(
        stale_target,
        column,
        R,
    ) == :stale_target_endpoint

    coefficient_tampering = merge(
        cert,
        (; witness = merge(witness, (; coefficient = 0))),
    )
    @test Suslin._validate_laurent_link_witness_certificate(
        coefficient_tampering,
        column,
        R,
    ) == :identity_replay_failed

    nonidentity = merge(
        cert,
        (; witness = merge(witness, (; exponent = (0, 0)))),
    )
    @test Suslin._validate_laurent_link_witness_certificate(
        nonidentity,
        column,
        R,
    ) == :identity_replay_failed
end

@testset "internal d14 Laurent link-witness candidate" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    source = _case008_d14_link_witness_source_data(fixture)
    report = case008_d14_laurent_link_witness_search_report(fixture)
    @test validate_case008_d14_laurent_link_witness_search_report(report, fixture) == :ok
    @test report.status == :candidate_found
    @test report.next_boundary == :laurent_link_witness_certificate

    candidate = first(report.candidates)
    witness = candidate.witness
    @test report.case_id == "case_008"
    @test witness.pivot_index == 10
    @test witness.partner_index == 1
    @test witness.coefficient == 1
    @test witness.exponent == (1, -1)

    @test Suslin._verify_laurent_link_witness_candidate(
        source.replay.after_column,
        fixture.ring,
        candidate,
    )
    cert = Suslin._laurent_link_witness_certificate_from_replay(
        source.context,
        witness,
        source.replay.after_column,
        fixture.ring,
    )
    @test cert.case_id == "case_008"
    @test cert.next_boundary == :laurent_endpoint_reduction
    @test Suslin._validate_laurent_link_witness_certificate(
        cert,
        source.replay.after_column,
        fixture.ring,
    ) == :ok
end
```

- [ ] **Step 2: Register the focused test**

In `test/runtests.jl`, add this entry immediately after
`"internal/laurent_descent_measure_helpers.jl",`:

```julia
        "internal/laurent_link_witness_helpers.jl",
```

- [ ] **Step 3: Verify the RED state**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_link_witness_helpers.jl")'
```

Expected: nonzero exit with an `UndefVarError` for one of the new
`Suslin._laurent_link_witness_*` helpers.

- [ ] **Step 4: Commit the failing contract**

Run:

```bash
git add test/internal/laurent_link_witness_helpers.jl test/runtests.jl
git commit -m "Add failing issue 339 internal helper contract"
```

---

### Task 2: Promote Generic Link-Witness Helpers

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Produces: `_LAURENT_LINK_WITNESS_FIELDS`
- Produces: `_LAURENT_LINK_WITNESS_CANDIDATE_FIELDS`
- Produces: `_laurent_link_witness_status(witness, n::Int, R)::Symbol`
- Produces: `_laurent_link_witness_operation(witness)::NamedTuple`
- Produces: `_laurent_link_endpoint_metadata(entry, R, entry_index::Int, column_measure; case_id)::NamedTuple`
- Produces: `_laurent_link_witness_candidate_from_replay(column, R, witness; case_id, require_strict = true)::NamedTuple`
- Produces: `_verify_laurent_link_witness_candidate(column, R, candidate)::Bool`

- [ ] **Step 1: Add constants after the Laurent descent certificate fields**

In `src/algorithm/column_reduction.jl`, after `_LAURENT_DESCENT_REQUIRED_MEASURE_FIELDS`, add:

```julia
const _LAURENT_LINK_WITNESS_FIELDS = (
    :family,
    :pivot_index,
    :partner_index,
    :coefficient,
    :exponent,
    :ring_generators,
)

const _LAURENT_LINK_WITNESS_CANDIDATE_FIELDS = (
    :witness,
    :source_endpoint,
    :target_endpoint,
    :replay_status,
    :identity_status,
    :measure_relation,
)
```

- [ ] **Step 2: Add status, operation, metadata, candidate, and verifier helpers**

After `_replay_laurent_elementary_entry_addition`, add:

```julia
function _laurent_link_witness_status(witness, n::Int, R)::Symbol
    _require_two_generator_laurent_ring(R)
    _laurent_descent_has_fields(witness, _LAURENT_LINK_WITNESS_FIELDS) ||
        return :malformed_witness
    witness.family == :two_entry_laurent_combination ||
        return :wrong_witness_family
    witness.ring_generators == _laurent_descent_ring_generators(R) ||
        return :wrong_ring_generators
    try
        pivot = _laurent_descent_checked_entry_index(
            witness.pivot_index,
            n,
            "pivot_index",
        )
        partner = _laurent_descent_checked_entry_index(
            witness.partner_index,
            n,
            "partner_index",
        )
        pivot != partner || return :pivot_partner_index_equality
        _laurent_descent_exponent_tuple(witness.exponent)
        _coerce_into_ring(R, witness.coefficient, "coefficient")
    catch err
        err isa InterruptException && rethrow()
        return :malformed_witness
    end
    return :ok
end

function _laurent_link_witness_operation(witness)
    return (;
        family = :entry_addition,
        target_index = witness.pivot_index,
        source_index = witness.partner_index,
        coefficient = witness.coefficient,
        exponent = witness.exponent,
        ring_generators = witness.ring_generators,
    )
end

function _laurent_link_endpoint_metadata(
    entry,
    R,
    entry_index::Int,
    column_measure;
    case_id,
)
    support = _laurent_descent_entry_support(_coerce_into_ring(R, entry, "entry"))
    return (;
        case_id,
        status = :link_witness_endpoint_metadata,
        entry_index,
        ring_generators = _laurent_descent_ring_generators(R),
        is_zero = isempty(support),
        term_count = length(support),
        support_bounds = _laurent_descent_support_bounds(support),
        leading_exponent = isempty(support) ? nothing : maximum(support),
        column_measure,
    )
end

function _laurent_link_witness_candidate_from_replay(
    column,
    R,
    witness;
    case_id,
    require_strict::Bool = true,
)
    status = _laurent_link_witness_status(witness, length(column), R)
    status == :ok ||
        throw(ArgumentError("invalid Laurent link witness: $(status)"))

    pivot = Int(witness.pivot_index)
    before_measure = _laurent_descent_measure_from_column(column, R; case_id)
    after_column = _replay_laurent_elementary_entry_addition(
        column,
        R,
        _laurent_link_witness_operation(witness),
    )
    after_measure = _laurent_descent_measure_from_column(after_column, R; case_id)
    relation = _strictly_decreases_laurent_measure(before_measure, after_measure)
    require_strict && !relation &&
        throw(ArgumentError("witness does not strictly decrease the replay measure"))
    return (;
        witness,
        source_endpoint = _laurent_link_endpoint_metadata(
            column[pivot],
            R,
            pivot,
            before_measure;
            case_id,
        ),
        target_endpoint = _laurent_link_endpoint_metadata(
            after_column[pivot],
            R,
            pivot,
            after_measure;
            case_id,
        ),
        replay_status = :ok,
        identity_status = :verified,
        measure_relation = relation ? :strict_decrease : :not_strict_decrease,
    )
end

function _verify_laurent_link_witness_candidate(column, R, candidate)::Bool
    try
        _laurent_descent_has_fields(
            candidate,
            _LAURENT_LINK_WITNESS_CANDIDATE_FIELDS,
        ) || return false
        candidate.replay_status == :ok || return false
        candidate.identity_status == :verified || return false
        candidate.measure_relation == :strict_decrease || return false

        witness = candidate.witness
        _laurent_link_witness_status(witness, length(column), R) == :ok ||
            return false
        expected = _laurent_link_witness_candidate_from_replay(
            column,
            R,
            witness;
            case_id = candidate.source_endpoint.case_id,
        )
        return candidate == expected
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 3: Verify GREEN for candidate helpers**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_link_witness_helpers.jl")'
```

Expected: the command still fails, but now at the missing certificate helper
stage rather than the candidate helper stage.

---

### Task 3: Promote Certificate Construction and Validation

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Produces: `_LAURENT_LINK_WITNESS_CERTIFICATE_STATUS`
- Produces: `_LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY`
- Produces: `_LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS`
- Produces: `_laurent_link_witness_certificate_from_replay(context, witness, column, R; require_strict = true)::NamedTuple`
- Produces: `_validate_laurent_link_witness_certificate(cert, column, R)::Symbol`

- [ ] **Step 1: Add certificate constants**

After `_LAURENT_LINK_WITNESS_CANDIDATE_FIELDS`, add:

```julia
const _LAURENT_LINK_WITNESS_CERTIFICATE_STATUS = :link_witness_certificate
const _LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY =
    :laurent_endpoint_reduction

const _LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS = (
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
```

- [ ] **Step 2: Add replay-derived certificate construction and validation**

After `_verify_laurent_link_witness_candidate`, add:

```julia
function _laurent_link_witness_certificate_from_replay(
    context,
    witness,
    column,
    R;
    require_strict::Bool = true,
)
    _laurent_descent_has_fields(
        context,
        (:case_id, :dimension, :ring_generators, :status),
    ) || throw(ArgumentError("context must include case_id, dimension, ring_generators, and status"))
    context.status == :link_witness_context ||
        throw(ArgumentError("context must have status :link_witness_context"))
    context.dimension == length(column) ||
        throw(ArgumentError("context dimension does not match source column"))
    context.ring_generators == _laurent_descent_ring_generators(R) ||
        throw(ArgumentError("context ring_generators do not match the ring"))

    candidate = _laurent_link_witness_candidate_from_replay(
        column,
        R,
        witness;
        case_id = context.case_id,
        require_strict,
    )
    _verify_laurent_link_witness_candidate(column, R, candidate) ||
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
        next_boundary = _LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY,
        status = _LAURENT_LINK_WITNESS_CERTIFICATE_STATUS,
    )
end

function _validate_laurent_link_witness_certificate(cert, column, R)::Symbol
    try
        _require_two_generator_laurent_ring(R)
        _laurent_descent_has_fields(
            cert,
            _LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS,
        ) || return :missing_certificate_fields
        cert.status == _LAURENT_LINK_WITNESS_CERTIFICATE_STATUS ||
            return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.identity_status == :verified || return :wrong_identity_status
        cert.next_boundary == _LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY ||
            return :wrong_next_boundary
        cert.dimension == length(column) || return :wrong_dimension
        cert.ring_generators == _laurent_descent_ring_generators(R) ||
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
        _verify_laurent_link_witness_candidate(column, R, candidate) ||
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

- [ ] **Step 3: Verify the internal focused command**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_link_witness_helpers.jl")'
```

Expected: exit 0.

- [ ] **Step 4: Commit the internal implementation**

Run:

```bash
git add src/algorithm/column_reduction.jl
git commit -m "Promote Laurent link witness helper internals"
```

---

### Task 4: Reuse Internal Helpers From Expert Tests

**Files:**
- Modify: `test/expert/case008_d14_laurent_link_witness_search.jl`
- Modify: `test/expert/laurent_link_witness_certificate.jl`

**Interfaces:**
- Preserves: `case008_d14_laurent_link_witness_search_report()`
- Preserves: `validate_case008_d14_laurent_link_witness_search_report(report, fixture)::Symbol`
- Preserves: `verify_laurent_link_witness_candidate(column, R, candidate)::Bool`
- Preserves: `laurent_link_witness_certificate(context, witness, column, R)::NamedTuple`
- Preserves: `validate_laurent_link_witness_certificate(cert, column, R)::Symbol`
- Preserves: `case008_d14_laurent_link_witness_certificate_summary()::NamedTuple`

- [ ] **Step 1: Replace generic search helpers with forwarding wrappers**

In `test/expert/case008_d14_laurent_link_witness_search.jl`, keep the d14
search constants and report fields. Replace the generic witness constants and
helper definitions with:

```julia
const LAURENT_LINK_WITNESS_FIELDS = Suslin._LAURENT_LINK_WITNESS_FIELDS
const LAURENT_LINK_WITNESS_CANDIDATE_FIELDS =
    Suslin._LAURENT_LINK_WITNESS_CANDIDATE_FIELDS

function _laurent_link_witness_has_fields(value, fields)::Bool
    return Suslin._laurent_descent_has_fields(value, fields)
end

function _laurent_link_witness_status(witness, n::Int, R)::Symbol
    return Suslin._laurent_link_witness_status(witness, n, R)
end

function _laurent_link_witness_operation(witness)
    return Suslin._laurent_link_witness_operation(witness)
end

function _laurent_link_witness_candidate_from_replay(
    column,
    R,
    witness;
    case_id,
    require_strict::Bool = true,
)
    return Suslin._laurent_link_witness_candidate_from_replay(
        column,
        R,
        witness;
        case_id,
        require_strict,
    )
end

function verify_laurent_link_witness_candidate(column, R, candidate)::Bool
    return Suslin._verify_laurent_link_witness_candidate(column, R, candidate)
end
```

- [ ] **Step 2: Replace certificate helpers with forwarding wrappers**

In `test/expert/laurent_link_witness_certificate.jl`, replace the local
certificate constants and generic helper implementations with:

```julia
const LAURENT_LINK_WITNESS_CERTIFICATE_STATUS =
    Suslin._LAURENT_LINK_WITNESS_CERTIFICATE_STATUS
const LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY =
    Suslin._LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY
const LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS =
    Suslin._LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS

function _laurent_link_witness_certificate_has_fields(value, fields)::Bool
    return Suslin._laurent_descent_has_fields(value, fields)
end

function _laurent_link_witness_certificate_from_replay(
    context,
    witness,
    column,
    R;
    require_strict::Bool = true,
)
    return Suslin._laurent_link_witness_certificate_from_replay(
        context,
        witness,
        column,
        R;
        require_strict,
    )
end

function laurent_link_witness_certificate(context, witness, column, R)
    return Suslin._laurent_link_witness_certificate_from_replay(
        context,
        witness,
        column,
        R;
        require_strict = true,
    )
end

function validate_laurent_link_witness_certificate(cert, column, R)::Symbol
    return Suslin._validate_laurent_link_witness_certificate(cert, column, R)
end
```

- [ ] **Step 3: Verify expert and internal focused commands**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_link_witness_helpers.jl")'
julia --project=. -e 'include("test/expert/laurent_link_witness_certificate.jl")'
```

Expected: both exit 0.

- [ ] **Step 4: Commit the expert reuse cleanup**

Run:

```bash
git add test/expert/case008_d14_laurent_link_witness_search.jl test/expert/laurent_link_witness_certificate.jl test/internal/laurent_link_witness_helpers.jl test/runtests.jl
git commit -m "Reuse internal Laurent link witness helpers"
```

---

### Task 5: Final Verification and PR Preparation

**Files:**
- No required edits.

**Interfaces:**
- Verifies all issue acceptance criteria and repository test gate.

- [ ] **Step 1: Run focused issue verification**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_link_witness_helpers.jl")'
julia --project=. -e 'include("test/expert/laurent_link_witness_certificate.jl")'
```

Expected: both exit 0.

- [ ] **Step 2: Run full required Agent Desk verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0. Existing Julia world-age warnings from Quillen fixture
catalog loading may appear, matching prior PRs.

- [ ] **Step 3: Run whitespace check**

Run:

```bash
git diff --check
```

Expected: exit 0.

- [ ] **Step 4: Review and create PR**

Run a whole-branch review, fix any Critical or Important findings, then use
`superpowers:finishing-a-development-branch`. Choose option 2, "Push and
create a Pull Request", under the Standing Answer Policy.
