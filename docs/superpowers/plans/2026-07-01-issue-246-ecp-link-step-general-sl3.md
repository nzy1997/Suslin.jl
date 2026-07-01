# Issue 246 ECP Link Step General SL3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Realize ECP link-step segments through verified ordinary-polynomial `SL_3` route certificates for #242 mainline link witnesses.

**Architecture:** Keep the existing fixture-backed link-step path as an explicit legacy mode, and add a general `:polynomial_sl3` route mode. General segments derive endpoint transport from verified endpoint reductions, route every elementary transport matrix through `_polynomial_factorization_route_certificate`, store the route certificates and factor groups, and verify route products plus endpoint maps.

**Tech Stack:** Julia, Oscar, existing Suslin ECP helpers in `src/algorithm/column_reduction.jl`, existing #184 route certificates in `src/algorithm/factorization.jl`, Julia `Test`.

## Global Constraints

- Preserve existing GF(2)/QQ fixture-backed link-step tests as regression coverage.
- At least one #242 mainline case must be realized without `:supplied_fixture_identity_sl2_endpoint_transport`, fixture ids, or GF(2)/QQ family recognizers as the support reason.
- Use `_polynomial_factorization_route_certificate` and `_verify_polynomial_factorization_route_certificate`; do not implement `SL_3` factorization in the ECP layer.
- The link-step verifier must accept non-identity `SL_2` blocks when backed by stored `SL_3` route certificates and must reject corrupted route factors, embedded blocks, endpoints, and factor ordering.
- General routed segments must store route matrices, route certificates, route factor groups, and metadata showing which embedded block was realized.
- Do not assemble the full recursive ECP reducer, implement #186 matrix peeling, or broaden Laurent/ToricBuilder support.
- Focused verification commands are `julia --project=. -e 'include("test/expert/ecp_link_step_general.jl")'` and `julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add `route_mode` storage, route-mode resolution, #184 route realization for endpoint elementary obligations, segment route metadata, and verifier checks.
- Create `test/expert/ecp_link_step_general.jl`: focused #246 #242 mainline route test and negative controls.
- Modify `test/expert/ecp_link_step.jl`: keep legacy fixture regression explicit with `route_mode = :legacy_fixture` where necessary.
- Modify `test/runtests.jl`: register the new expert test.
- Add this plan and the committed design spec.

### Task 1: Add General Link-Step Route Tests

**Files:**
- Create: `test/expert/ecp_link_step_general.jl`
- Modify: `test/runtests.jl`
- Modify: `test/expert/ecp_link_step.jl`

**Interfaces:**
- Consumes future `ecp_link_step_certificate(...; route_mode = :polynomial_sl3)`.
- Produces focused RED coverage for route-backed multi-segment link steps, route certificate verification, endpoint mapping, and tamper rejection.

- [ ] **Step 1: Create the failing focused test**

Create `test/expert/ecp_link_step_general.jl`:

```julia
using Test
using Oscar
using Suslin

include(joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl"))

function _general_link_column(entry)
    return [getproperty(entry.column_entries, name) for name in entry.column_order]
end

function _general_link_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _general_link_apply(factors, column, R)
    return _general_link_factor_product(factors, R, length(column)) *
           matrix(R, length(column), 1, collect(column))
end

function _general_replace_field(record, constructor, field::Symbol, value)
    fields = fieldnames(typeof(record))
    values = [getfield(record, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown record field: $(field)")
    values[idx] = value
    return constructor(values...)
end

function _general_replace_segment_field(segment::NamedTuple, field::Symbol, value)
    haskey(segment, field) || error("unknown segment field: $(field)")
    return merge(segment, NamedTuple{(field,)}((value,)))
end

function _general_tamper_segment(record, segment_index::Int, field::Symbol, value)
    segments = collect(record.segments)
    segments[segment_index] = _general_replace_segment_field(segments[segment_index], field, value)
    return _general_replace_field(record, Suslin.ECPLinkStepCertificate, :segments, tuple(segments...))
end

function _general_replace_route_certificate_factor(certificate, factor_index::Int, replacement)
    factors = copy(certificate.factors)
    factors[factor_index] = replacement
    return Suslin.PolynomialFactorizationRouteCertificate(
        certificate.matrix,
        certificate.route,
        factors,
        certificate.product,
        certificate.evidence,
        certificate.status,
        certificate.verification,
    )
end

@testset "ECP link steps route through general SL3 certificates" begin
    cases = ECPMainlineFixtureCatalog.cases_by_id()
    entry = cases["ecp-mainline-sl3-route-qq"]
    R = entry.ring.object
    column = _general_link_column(entry)
    witness = entry.support_evidence.link_witness

    record = Suslin.ecp_link_step_certificate(
        column,
        R;
        link_witness = witness,
        route_mode = :polynomial_sl3,
    )

    @test record isa Suslin.ECPLinkStepCertificate
    @test record.route_mode == :polynomial_sl3
    @test Suslin.verify_ecp_link_step_certificate(record)
    @test length(record.segments) >= 2
    @test all(segment -> segment.support_family == :polynomial_sl3_route_endpoint_transport, record.segments)
    @test all(segment -> segment.support_family != :supplied_fixture_identity_sl2_endpoint_transport, record.segments)
    @test any(segment -> segment.sl2_block != identity_matrix(R, 2), record.segments)

    for segment in record.segments
        @test !isempty(segment.sl3_route_certificates)
        @test length(segment.sl3_route_certificates) == length(segment.sl3_route_matrices)
        @test length(segment.sl3_route_factor_groups) == length(segment.sl3_route_certificates)
        @test length(segment.sl3_route_metadata) == length(segment.sl3_route_certificates)
        for idx in eachindex(segment.sl3_route_certificates)
            route_cert = segment.sl3_route_certificates[idx]
            @test route_cert isa Suslin.PolynomialFactorizationRouteCertificate
            @test route_cert.status == :supported
            @test route_cert.matrix == segment.sl3_route_matrices[idx]
            @test Suslin._verify_polynomial_factorization_route_certificate(route_cert)
            @test segment.sl3_route_factor_groups[idx] == tuple(route_cert.factors...)
            @test segment.sl3_route_metadata[idx].route == route_cert.route
            @test segment.sl3_route_metadata[idx].factor_count == length(route_cert.factors)
        end
        @test _general_link_apply(segment.forward_factors, segment.from_column, R) ==
              matrix(R, length(segment.to_column), 1, collect(segment.to_column))
        @test _general_link_apply(segment.inverse_factors, segment.to_column, R) ==
              matrix(R, length(segment.from_column), 1, collect(segment.from_column))
    end

    @test _general_link_apply(record.forward_factors, record.lower_variable_column, R) ==
          matrix(R, length(record.transformed_column), 1, collect(record.transformed_column))
    @test _general_link_apply(record.reduction_factors, record.transformed_column, R) ==
          matrix(R, length(record.lower_variable_column), 1, collect(record.lower_variable_column))
end
```

- [ ] **Step 2: Add negative controls**

Append to the same testset after the positive assertions:

```julia
    first_segment = record.segments[1]
    first_route = first_segment.sl3_route_certificates[1]
    corrupted_route = _general_replace_route_certificate_factor(
        first_route,
        1,
        elementary_matrix(nrows(first_route.factors[1]), 1, 2, one(R), R),
    )
    corrupted_routes = collect(first_segment.sl3_route_certificates)
    corrupted_routes[1] = corrupted_route
    @test !Suslin.verify_ecp_link_step_certificate(
        _general_tamper_segment(record, 1, :sl3_route_certificates, tuple(corrupted_routes...)),
    )

    bad_sl2 = first_segment.sl2_block == identity_matrix(R, 2) ?
        elementary_matrix(2, 1, 2, one(R), R) :
        identity_matrix(R, 2)
    @test !Suslin.verify_ecp_link_step_certificate(_general_tamper_segment(record, 1, :sl2_block, bad_sl2))

    bad_endpoint = ntuple(
        idx -> idx == 1 ? first_segment.to_column[idx] + one(R) : first_segment.to_column[idx],
        length(first_segment.to_column),
    )
    @test !Suslin.verify_ecp_link_step_certificate(_general_tamper_segment(record, 1, :to_column, bad_endpoint))

    if length(record.forward_factors) >= 2
        reordered = copy(record.forward_factors)
        reordered[1], reordered[2] = reordered[2], reordered[1]
        @test !Suslin.verify_ecp_link_step_certificate(
            _general_replace_field(record, Suslin.ECPLinkStepCertificate, :forward_factors, reordered),
        )
    end
```

- [ ] **Step 3: Register the test**

In `test/runtests.jl`, add:

```julia
"expert/ecp_link_step_general.jl",
```

next to `"expert/ecp_link_step.jl"`.

- [ ] **Step 4: Keep fixture regression explicit**

In `test/expert/ecp_link_step.jl`, update the GF(2) and QQ positive
`Suslin.ecp_link_step_certificate` calls to include:

```julia
route_mode = :legacy_fixture,
```

so those assertions continue to cover the legacy family tag deliberately.

- [ ] **Step 5: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step_general.jl")'
```

Expected: FAIL with `UndefKeywordError: keyword argument route_mode not assigned` or equivalent missing-field failures, because the production route mode is not implemented yet.

### Task 2: Add Route-Backed Segment Realization

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/ecp_link_step_general.jl`
- Test: `test/expert/ecp_link_step.jl`

**Interfaces:**
- Consumes `_polynomial_factorization_route_certificate`, `_verify_polynomial_factorization_route_certificate`, existing endpoint reductions, and existing link identity replay.
- Produces `route_mode` storage, route-backed segment construction, and segment metadata fields `sl3_route_matrices`, `sl3_route_certificates`, `sl3_route_factor_groups`, and `sl3_route_metadata`.

- [ ] **Step 1: Add route mode to the certificate and constructor**

In `ECPLinkStepCertificate`, insert `route_mode::Symbol` after `path_columns`.
In `ecp_link_step_certificate`, add the keyword:

```julia
route_mode::Symbol = :auto,
```

Resolve it before segment construction:

```julia
resolved_route_mode = _ecp_link_step_resolve_route_mode(witness, route_mode)
segments = _ecp_link_step_segments(witness, path_columns; route_mode = resolved_route_mode)
```

Pass `resolved_route_mode` into both `ECPLinkStepCertificate` constructors.

- [ ] **Step 2: Add route-mode helpers**

Add near the existing link-step helpers:

```julia
function _ecp_link_step_resolve_route_mode(witness::ECPLinkWitnessRecord, route_mode::Symbol)
    route_mode in (:auto, :legacy_fixture, :polynomial_sl3) ||
        throw(ArgumentError("unsupported ECP link step route_mode $(route_mode)"))
    if route_mode == :auto
        return (_ecp_link_step_matches_gf2_fixture(witness) ||
                _ecp_link_step_matches_qq_fixture(witness)) ?
            :legacy_fixture :
            :polynomial_sl3
    end
    return route_mode
end

function _ecp_link_step_supported_family(witness::ECPLinkWitnessRecord, route_mode::Symbol)
    resolved = _ecp_link_step_resolve_route_mode(witness, route_mode)
    if resolved == :legacy_fixture
        if _ecp_link_step_matches_gf2_fixture(witness) || _ecp_link_step_matches_qq_fixture(witness)
            return :supplied_fixture_identity_sl2_endpoint_transport
        end
        probe_ids = tuple((probe.id for probe in witness.residue_probes)...)
        throw(ArgumentError("unsupported ECP legacy fixture link step family for supplied link witness probes $(probe_ids)"))
    end
    return :polynomial_sl3_route_endpoint_transport
end
```

Change `_ecp_link_step_segments` and `_ecp_link_step_segment` signatures so
they accept `; route_mode::Symbol = :auto` and use
`_ecp_link_step_supported_family(witness, route_mode)`.

- [ ] **Step 3: Route endpoint elementary obligations through #184**

Replace the old transport call with `_ecp_link_step_endpoint_transport`. The
general helper should compute endpoint reductions as before, then route each
raw elementary matrix:

```julia
function _ecp_link_step_endpoint_transport(R, from_column, to_column, link_identity, support_family::Symbol)
    link_identity.overall_ok ||
        throw(ArgumentError("ECP link step requires a replayed link identity"))
    from_certificate = try
        ecp_column_reduction_certificate(collect(from_column), R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("unsupported ECP link step source path column"))
    end
    to_certificate = try
        ecp_column_reduction_certificate(collect(to_column), R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("unsupported ECP link step target path column"))
    end
    raw_factors = vcat(_ecp_inverse_factor_sequence(to_certificate.factors), from_certificate.factors)

    if support_family == :supplied_fixture_identity_sl2_endpoint_transport
        endpoint_transport_matrix = _factor_sequence_product(raw_factors, R, length(from_column))
        return (; from_certificate, to_certificate, endpoint_transport_matrix, factors = raw_factors,
            sl3_route_matrices = (), sl3_route_certificates = (), sl3_route_factor_groups = (),
            sl3_route_metadata = ())
    elseif support_family == :polynomial_sl3_route_endpoint_transport
        sl3_route_matrices = tuple(raw_factors...)
        sl3_route_certificates = tuple((
            _ecp_link_step_route_certificate(factor)
            for factor in sl3_route_matrices
        )...)
        sl3_route_factor_groups = tuple((tuple(cert.factors...) for cert in sl3_route_certificates)...)
        factors = Any[]
        for group in sl3_route_factor_groups
            append!(factors, group)
        end
        endpoint_transport_matrix = _factor_sequence_product(factors, R, length(from_column))
        sl3_route_metadata = _ecp_link_step_route_metadata(sl3_route_certificates)
        return (; from_certificate, to_certificate, endpoint_transport_matrix, factors,
            sl3_route_matrices, sl3_route_certificates, sl3_route_factor_groups,
            sl3_route_metadata)
    end
    throw(ArgumentError("unsupported ECP link step family $(support_family)"))
end

function _ecp_link_step_route_certificate(factor)
    certificate = _polynomial_factorization_route_certificate(factor; allow_recursive_column_peel = false)
    certificate.status == :supported ||
        throw(ArgumentError("ECP link step SL_3 route obligation is staged"))
    _verify_polynomial_factorization_route_certificate(certificate) ||
        throw(ArgumentError("ECP link step SL_3 route certificate does not verify"))
    certificate.product == factor ||
        throw(ArgumentError("ECP link step SL_3 route certificate product does not match its obligation"))
    return certificate
end

function _ecp_link_step_route_metadata(route_certificates)
    return tuple((
        (;
            source = :polynomial_factorization_route_certificate,
            obligation_index = idx,
            route = route_certificates[idx].route,
            status = route_certificates[idx].status,
            factor_count = length(route_certificates[idx].factors),
        )
        for idx in eachindex(route_certificates)
    )...)
end
```

After computing `transport`, set:

```julia
sl2_block = _ecp_link_step_embedded_sl2_block(
    transport.sl3_route_matrices,
    R,
    n,
    (1, 2),
)
sl2_embedding = block_embedding(sl2_block, n, (1, 2))
```

and add all `transport.sl3_*` fields to the returned segment.

- [ ] **Step 4: Extract the visible embedded SL2 block**

Add:

```julia
function _ecp_link_step_embedded_sl2_block(route_matrices, R, n::Int, indices)
    identity_block = identity_matrix(R, 2)
    for route_matrix in route_matrices
        block = _ecp_link_step_extract_embedded_sl2_block(route_matrix, R, n, indices)
        block === nothing && continue
        block == identity_block && continue
        return block
    end
    return identity_block
end

function _ecp_link_step_extract_embedded_sl2_block(route_matrix, R, n::Int, indices)
    nrows(route_matrix) == n && ncols(route_matrix) == n || return nothing
    _same_base_ring(base_ring(route_matrix), R) || return nothing
    i, j = indices
    block = matrix(R, 2, 2, [route_matrix[i, i], route_matrix[i, j], route_matrix[j, i], route_matrix[j, j]])
    return block_embedding(block, n, indices) == route_matrix ? block : nothing
end
```

- [ ] **Step 5: Run GREEN focused tests**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step_general.jl")'
julia --project=. -e 'include("test/expert/ecp_link_step.jl")'
```

Expected: both PASS. If failures occur, use `superpowers:systematic-debugging`
before changing production code again.

### Task 3: Harden Replay Verification and Final Integration

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/runtests.jl`
- Test: `test/expert/ecp_link_step_general.jl`
- Test: `test/expert/park_woodburn_route_certificate.jl`

**Interfaces:**
- Consumes Task 2 segment fields.
- Produces verifier checks that route certificates, route products, route factor groups, metadata, non-identity embedded blocks, endpoint maps, composed factors, and stored verification tuples all reject tampering.

- [ ] **Step 1: Verify route fields inside `_ecp_link_step_segment_verification`**

Extend the function signature with:

```julia
sl3_route_matrices,
sl3_route_certificates,
sl3_route_factor_groups,
sl3_route_metadata,
```

Compute:

```julia
legacy_support = support_family == :supplied_fixture_identity_sl2_endpoint_transport
general_support = support_family == :polynomial_sl3_route_endpoint_transport
sl2_identity_ok = sl2_block == identity_matrix(R, 2)
sl2_block_ok = nrows(sl2_block) == 2 &&
    ncols(sl2_block) == 2 &&
    _same_base_ring(base_ring(sl2_block), R) &&
    det(sl2_block) == one(R) &&
    (legacy_support ? sl2_identity_ok : true)
sl2_embedding_ok = sl2_embedding == block_embedding(sl2_block, n, (1, 2))
route_certificates_ok = _ecp_link_step_route_certificates_ok(
    sl3_route_matrices,
    sl3_route_certificates,
    sl3_route_factor_groups,
    sl3_route_metadata,
    R,
    n,
    support_family,
)
embedded_sl2_route_ok = legacy_support ||
    _ecp_link_step_has_embedded_sl2_route(sl3_route_matrices, sl2_embedding)
```

Include `route_certificates_ok` and `embedded_sl2_route_ok` in `overall_ok` and
return them in the verification tuple.

- [ ] **Step 2: Add route verification helpers**

Add:

```julia
function _ecp_link_step_route_certificates_ok(
    sl3_route_matrices,
    sl3_route_certificates,
    sl3_route_factor_groups,
    sl3_route_metadata,
    R,
    n::Int,
    support_family::Symbol,
)::Bool
    if support_family == :supplied_fixture_identity_sl2_endpoint_transport
        return isempty(sl3_route_matrices) &&
            isempty(sl3_route_certificates) &&
            isempty(sl3_route_factor_groups) &&
            isempty(sl3_route_metadata)
    end
    support_family == :polynomial_sl3_route_endpoint_transport || return false
    length(sl3_route_matrices) == length(sl3_route_certificates) || return false
    length(sl3_route_factor_groups) == length(sl3_route_certificates) || return false
    length(sl3_route_metadata) == length(sl3_route_certificates) || return false
    isempty(sl3_route_certificates) && return false
    for idx in eachindex(sl3_route_certificates)
        matrix_i = sl3_route_matrices[idx]
        cert = sl3_route_certificates[idx]
        nrows(matrix_i) == n && ncols(matrix_i) == n || return false
        _same_base_ring(base_ring(matrix_i), R) || return false
        cert isa PolynomialFactorizationRouteCertificate || return false
        cert.status == :supported || return false
        cert.matrix == matrix_i || return false
        cert.product == matrix_i || return false
        _verify_polynomial_factorization_route_certificate(cert) || return false
        tuple(cert.factors...) == sl3_route_factor_groups[idx] || return false
        metadata = sl3_route_metadata[idx]
        hasproperty(metadata, :source) && metadata.source == :polynomial_factorization_route_certificate || return false
        hasproperty(metadata, :obligation_index) && metadata.obligation_index == idx || return false
        hasproperty(metadata, :route) && metadata.route == cert.route || return false
        hasproperty(metadata, :status) && metadata.status == cert.status || return false
        hasproperty(metadata, :factor_count) && metadata.factor_count == length(cert.factors) || return false
    end
    return true
end

function _ecp_link_step_has_embedded_sl2_route(sl3_route_matrices, sl2_embedding)::Bool
    return any(route_matrix -> route_matrix == sl2_embedding, sl3_route_matrices)
end
```

- [ ] **Step 3: Include route fields in segment equivalence**

In `_ecp_link_step_segment_equivalent`, add equality checks for:

```julia
left.sl3_route_matrices == right.sl3_route_matrices &&
left.sl3_route_certificates == right.sl3_route_certificates &&
left.sl3_route_factor_groups == right.sl3_route_factor_groups &&
left.sl3_route_metadata == right.sl3_route_metadata &&
```

In `_ecp_link_step_replay_summary`, recompute with:

```julia
_ecp_link_step_segments(certificate.link_witness, recomputed_path_columns; route_mode = certificate.route_mode)
```

and require `certificate.route_mode` to be `:legacy_fixture` or
`:polynomial_sl3`.

- [ ] **Step 4: Run required focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step_general.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: both PASS.

- [ ] **Step 5: Run full package verification and commit**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: both PASS. Commit with:

```bash
git add src/algorithm/column_reduction.jl test/expert/ecp_link_step_general.jl test/expert/ecp_link_step.jl test/runtests.jl
git commit -m "feat: route ecp link steps through sl3 certificates"
```

## Plan Self-Review

- Every issue acceptance requirement maps to at least one task and focused assertion.
- The route API is reused; no local `SL_3` factorization is specified.
- The legacy fixture path remains explicit and testable.
- Unsupported route obligations fail with `ArgumentError` rather than guessed data.
- No incomplete markers or contradictory task steps remain.
