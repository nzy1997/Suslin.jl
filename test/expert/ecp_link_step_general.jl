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

function _general_link_witness_data(witness; residue_probes = witness.residue_probes)
    return (;
        source = :supplied_link_witness,
        residue_probes,
        tail_reductions = witness.tail_reductions,
        resultants = witness.resultants,
        bezout_coefficients = witness.bezout_coefficients,
        coverage_multipliers = witness.coverage_multipliers,
        path_points = witness.path_points,
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
        matching_route_metadata = nothing
        for idx in eachindex(segment.sl3_route_certificates)
            route_cert = segment.sl3_route_certificates[idx]
            metadata = segment.sl3_route_metadata[idx]
            @test route_cert isa Suslin.PolynomialFactorizationRouteCertificate
            @test route_cert.status == :supported
            @test route_cert.matrix == segment.sl3_route_matrices[idx]
            @test Suslin._verify_polynomial_factorization_route_certificate(route_cert)
            @test segment.sl3_route_factor_groups[idx] == tuple(route_cert.factors...)
            @test metadata.route == route_cert.route
            @test metadata.factor_count == length(route_cert.factors)
            if route_cert.matrix == segment.sl2_embedding
                matching_route_metadata = metadata
                @test metadata.embedded_block_indices == (1, 2)
                @test metadata.embedded_block_matrix == segment.sl2_block
            end
        end
        @test !isnothing(matching_route_metadata)
        @test _general_link_apply(segment.forward_factors, segment.from_column, R) ==
              matrix(R, length(segment.to_column), 1, collect(segment.to_column))
        @test _general_link_apply(segment.inverse_factors, segment.to_column, R) ==
              matrix(R, length(segment.from_column), 1, collect(segment.from_column))
    end

    @test _general_link_apply(record.forward_factors, record.lower_variable_column, R) ==
          matrix(R, length(record.transformed_column), 1, collect(record.transformed_column))
    @test _general_link_apply(record.reduction_factors, record.transformed_column, R) ==
          matrix(R, length(record.lower_variable_column), 1, collect(record.lower_variable_column))

    nonfixture_probes = (
        merge(witness.residue_probes[1], (; maximal_ideal_generators = (entry.ring.generators[1],))),
        witness.residue_probes[2],
    )
    nonfixture_witness = Suslin.ecp_link_witness(
        column,
        R;
        variable_order = entry.ring.generators,
        selected_variable = entry.selected_variable.generator,
        supplied_link_witness = _general_link_witness_data(
            witness;
            residue_probes = nonfixture_probes,
        ),
    )
    @test Suslin.verify_ecp_link_witness(nonfixture_witness)
    @test !Suslin._ecp_link_step_matches_qq_fixture(nonfixture_witness)
    auto_record = Suslin.ecp_link_step_certificate(column, R; link_witness = nonfixture_witness)
    @test auto_record.route_mode == :polynomial_sl3
    @test all(
        segment -> segment.support_family == :polynomial_sl3_route_endpoint_transport,
        auto_record.segments,
    )
    @test Suslin.verify_ecp_link_step_certificate(auto_record)

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
end
