using Test
using Oscar
using Suslin

const ECP_COLUMN_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl")
include(ECP_COLUMN_CATALOG_PATH)

function _case_by_id(id::AbstractString)
    return ECPColumnFixtureCatalog.cases_by_id()[id]
end

function _column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _gf2_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    v = _column(entry)
    G = y * v[2] + v[3]
    return (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((;
            probe_id = :gf2_fixture_probe,
            G,
            lifted_tail_coefficients = (y, one(R)),
            tilde_G = G,
        ),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
end

function _qq_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    return (;
        source = :supplied_link_witness,
        residue_probes = (
            (; id = :qq_y_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),
            (; id = :qq_x_probe, kind = :deterministic_fixture, maximal_ideal_generators = (x,)),
        ),
        tail_reductions = (
            (; probe_id = :qq_y_probe, G = y, lifted_tail_coefficients = (zero(R), one(R)), tilde_G = y),
            (; probe_id = :qq_x_probe, G = x, lifted_tail_coefficients = (one(R), zero(R)), tilde_G = x),
        ),
        resultants = (y^2, y + one(R)),
        bezout_coefficients = (
            (; f = zero(R), h = y),
            (; f = one(R), h = -x),
        ),
        coverage_multipliers = (one(R), one(R) - y),
        path_points = (zero(R), y^2 * x, x),
    )
end

function _mutate_witness(
    witness;
    residue_probes = witness.residue_probes,
    tail_reductions = witness.tail_reductions,
    resultants = witness.resultants,
    bezout_coefficients = witness.bezout_coefficients,
    coverage_multipliers = witness.coverage_multipliers,
    path_points = witness.path_points,
)
    return merge(witness, (;
        residue_probes,
        tail_reductions,
        resultants,
        bezout_coefficients,
        coverage_multipliers,
        path_points,
    ))
end

function _replace_record_field(record, constructor, field::Symbol, value)
    fields = fieldnames(typeof(record))
    values = [getfield(record, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown record field: $(field)")
    values[idx] = value
    return constructor(values...)
end

function _replace_segment_field(segment::NamedTuple, field::Symbol, value)
    haskey(segment, field) || error("unknown segment field: $(field)")
    return merge(segment, NamedTuple{(field,)}((value,)))
end

function _replace_segment_field(segment, field::Symbol, value)
    fields = fieldnames(typeof(segment))
    values = [getfield(segment, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown segment field: $(field)")
    values[idx] = value
    return typeof(segment)(values...)
end

function _tamper_segment(record, segment_index::Int, field::Symbol, value)
    segments = collect(record.segments)
    segments[segment_index] = _replace_segment_field(segments[segment_index], field, value)
    return _replace_record_field(record, Suslin.ECPLinkStepCertificate, :segments, tuple(segments...))
end

function _link_step_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _link_step_apply_factors(factors, column, R)
    return _link_step_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _assert_link_step_segment_maps(record)
    R = record.ring
    for segment in record.segments
        @test _link_step_apply_factors(segment.forward_factors, segment.from_column, R) ==
              matrix(R, length(segment.to_column), 1, collect(segment.to_column))
        @test _link_step_apply_factors(segment.inverse_factors, segment.to_column, R) ==
              matrix(R, length(segment.from_column), 1, collect(segment.from_column))
    end
end

function _assert_link_step_composed_maps(record)
    R = record.ring
    @test _link_step_apply_factors(record.forward_factors, record.lower_variable_column, R) ==
          matrix(R, length(record.transformed_column), 1, collect(record.transformed_column))
    @test _link_step_apply_factors(record.reduction_factors, record.transformed_column, R) ==
          matrix(R, length(record.lower_variable_column), 1, collect(record.lower_variable_column))
end

function _unsupported_identity_sl2_case()
    R, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    column = [x^2 + y + one(R), x * y^2, x^2 + y^2]
    r1 = resultant(column[1], column[2], 1)
    r2 = resultant(column[1], column[3], 1)
    bezout1 = coordinates(r1, ideal(R, [column[1], column[2]]))
    bezout2 = coordinates(r2, ideal(R, [column[1], column[3]]))
    coverage = coordinates(one(R), ideal(R, [r1, r2]))
    g1 = coverage[1, 1]
    g2 = coverage[1, 2]
    witness = (;
        source = :supplied_link_witness,
        residue_probes = (
            (; id = :unsupported_p1, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),
            (; id = :unsupported_p2, kind = :deterministic_fixture, maximal_ideal_generators = (x,)),
        ),
        tail_reductions = (
            (; probe_id = :unsupported_p1, G = column[2], lifted_tail_coefficients = (one(R), zero(R)), tilde_G = column[2]),
            (; probe_id = :unsupported_p2, G = column[3], lifted_tail_coefficients = (zero(R), one(R)), tilde_G = column[3]),
        ),
        resultants = (r1, r2),
        bezout_coefficients = (
            (; f = bezout1[1, 1], h = bezout1[1, 2]),
            (; f = bezout2[1, 1], h = bezout2[1, 2]),
        ),
        coverage_multipliers = (g1, g2),
        path_points = (zero(R), r1 * g1 * x, x),
    )
    return (; R, x, y, column, witness)
end

function _relabelled_unit_tail_case()
    R, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    column = [x + y, one(R), zero(R)]
    witness = (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((;
            probe_id = :gf2_fixture_probe,
            G = one(R),
            lifted_tail_coefficients = (one(R), zero(R)),
            tilde_G = one(R),
        ),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = zero(R), h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
    return (; R, x, y, column, witness)
end

function _extra_generator_gf2_case()
    R, (x, y, z) = Oscar.polynomial_ring(GF(2), ["x", "y", "z"])
    column = [x + y^2, x * y + x + one(R), x^2 + x * y + y + one(R)]
    G = y * column[2] + column[3]
    witness = (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((;
            probe_id = :gf2_fixture_probe,
            G,
            lifted_tail_coefficients = (y, one(R)),
            tilde_G = G,
        ),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
    return (; R, x, y, z, column, witness)
end

@testset "ECP link step certificate replays path transport" begin
    gf2_entry = _case_by_id("ecp-variable-change-monic-gf2")
    gf2_column = _column(gf2_entry)
    gf2_record = Suslin.ecp_link_step_certificate(
        gf2_column,
        gf2_entry.ring.object;
        variable_order = gf2_entry.ring.generators,
        selected_variable = gf2_entry.ring.generators[1],
        supplied_link_witness = _gf2_link_witness(gf2_entry),
    )
    @test gf2_record.verification.overall_ok == true
    @test length(gf2_record.path_columns) == 2
    @test length(gf2_record.segments) == 1
    @test gf2_record.lower_variable_column == gf2_record.path_columns[1]
    @test gf2_record.transformed_column == tuple(gf2_column...)
    @test all(segment -> segment.support_family == :supplied_fixture_identity_sl2_endpoint_transport, gf2_record.segments)
    @test all(segment -> segment.verification.endpoint_transport_ok, gf2_record.segments)
    _assert_link_step_segment_maps(gf2_record)
    _assert_link_step_composed_maps(gf2_record)
    @test Suslin.verify_ecp_link_step_certificate(gf2_record)

    qq_entry = _case_by_id("ecp-monic-first-entry-qq")
    qq_column = _column(qq_entry)
    qq_witness = _qq_link_witness(qq_entry)
    qq_record = Suslin.ecp_link_step_certificate(
        qq_column,
        qq_entry.ring.object;
        variable_order = qq_entry.ring.generators,
        selected_variable = qq_entry.ring.generators[1],
        supplied_link_witness = qq_witness,
    )
    @test qq_record.verification.overall_ok == true
    @test length(qq_record.path_columns) == 3
    @test length(qq_record.segments) == 2
    @test qq_record.lower_variable_column == qq_record.path_columns[1]
    @test qq_record.transformed_column == tuple(qq_column...)
    @test any(segment -> segment.delta != zero(qq_entry.ring.object), qq_record.segments)
    @test qq_record.verification.composed_forward_ok == true
    @test qq_record.verification.composed_reduction_ok == true
    @test all(segment -> segment.support_family == :supplied_fixture_identity_sl2_endpoint_transport, qq_record.segments)
    @test all(segment -> segment.verification.endpoint_transport_ok, qq_record.segments)
    _assert_link_step_segment_maps(qq_record)
    _assert_link_step_composed_maps(qq_record)
    @test Suslin.verify_ecp_link_step_certificate(qq_record)

    R = qq_entry.ring.object
    x, y = qq_entry.ring.generators
    bad_sl2 = identity_matrix(R, 2) + elementary_matrix(2, 1, 2, one(R), R)
    @test !Suslin.verify_ecp_link_step_certificate(_tamper_segment(qq_record, 1, :sl2_block, bad_sl2))
    @test !Suslin.verify_ecp_link_step_certificate(_tamper_segment(
        qq_record,
        1,
        :elementary_factors,
        qq_record.segments[1].elementary_factors[1:(end - 1)],
    ))
    @test_throws ArgumentError Suslin.ecp_link_step_certificate(
        qq_column,
        R;
        variable_order = qq_entry.ring.generators,
        selected_variable = x,
        supplied_link_witness = _mutate_witness(
            qq_witness;
            path_points = (zero(R), x * y^2 + one(R), x),
        ),
    )

    unsupported = _unsupported_identity_sl2_case()
    unsupported_witness = Suslin.ecp_link_witness(
        unsupported.column,
        unsupported.R;
        variable_order = (unsupported.x, unsupported.y),
        selected_variable = unsupported.x,
        supplied_link_witness = unsupported.witness,
    )
    @test Suslin.verify_ecp_link_witness(unsupported_witness)
    @test_throws ArgumentError Suslin.ecp_link_step_certificate(
        unsupported.column,
        unsupported.R;
        link_witness = unsupported_witness,
    )

    relabeled_witness = _mutate_witness(
        unsupported.witness;
        residue_probes = (
            merge(unsupported.witness.residue_probes[1], (; id = :qq_y_probe)),
            merge(unsupported.witness.residue_probes[2], (; id = :qq_x_probe)),
        ),
        tail_reductions = (
            merge(unsupported.witness.tail_reductions[1], (; probe_id = :qq_y_probe)),
            merge(unsupported.witness.tail_reductions[2], (; probe_id = :qq_x_probe)),
        ),
    )
    relabeled_record = Suslin.ecp_link_witness(
        unsupported.column,
        unsupported.R;
        variable_order = (unsupported.x, unsupported.y),
        selected_variable = unsupported.x,
        supplied_link_witness = relabeled_witness,
    )
    @test Suslin.verify_ecp_link_witness(relabeled_record)
    @test_throws ArgumentError Suslin.ecp_link_step_certificate(
        unsupported.column,
        unsupported.R;
        link_witness = relabeled_record,
    )

    relabelled_unit_tail = _relabelled_unit_tail_case()
    relabelled_unit_tail_record = Suslin.ecp_link_witness(
        relabelled_unit_tail.column,
        relabelled_unit_tail.R;
        variable_order = (relabelled_unit_tail.x, relabelled_unit_tail.y),
        selected_variable = relabelled_unit_tail.x,
        supplied_link_witness = relabelled_unit_tail.witness,
    )
    @test Suslin.verify_ecp_link_witness(relabelled_unit_tail_record)
    @test_throws ArgumentError Suslin.ecp_link_step_certificate(
        relabelled_unit_tail.column,
        relabelled_unit_tail.R;
        link_witness = relabelled_unit_tail_record,
    )

    extra_generator = _extra_generator_gf2_case()
    extra_generator_record = Suslin.ecp_link_witness(
        extra_generator.column,
        extra_generator.R;
        variable_order = (extra_generator.x, extra_generator.y, extra_generator.z),
        selected_variable = extra_generator.x,
        supplied_link_witness = extra_generator.witness,
    )
    @test Suslin.verify_ecp_link_witness(extra_generator_record)
    @test_throws ArgumentError Suslin.ecp_link_step_certificate(
        extra_generator.column,
        extra_generator.R;
        link_witness = extra_generator_record,
    )
end
