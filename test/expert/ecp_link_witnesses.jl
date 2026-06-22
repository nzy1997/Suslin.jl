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

function _replace_tuple_entry(values::Tuple, idx::Int, value)
    return ntuple(j -> j == idx ? value : values[j], length(values))
end

function _mutate_witness(
    witness;
    resultants = witness.resultants,
    bezout_coefficients = witness.bezout_coefficients,
    coverage_multipliers = witness.coverage_multipliers,
    path_points = witness.path_points,
)
    return merge(witness, (;
        resultants,
        bezout_coefficients,
        coverage_multipliers,
        path_points,
    ))
end

function _tamper_record_field(record, field::Symbol, value)
    fields = fieldnames(typeof(record))
    values = [getfield(record, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown record field: $(field)")
    values[idx] = value
    return Suslin.ECPLinkWitnessRecord(values...)
end

@testset "supplied ECP link witnesses replay exactly" begin
    gf2_entry = _case_by_id("ecp-variable-change-monic-gf2")
    gf2_column = _column(gf2_entry)
    gf2_record = Suslin.ecp_link_witness(
        gf2_column,
        gf2_entry.ring.object;
        variable_order = gf2_entry.ring.generators,
        selected_variable = gf2_entry.ring.generators[1],
        supplied_link_witness = _gf2_link_witness(gf2_entry),
    )
    @test !any(is_unit, gf2_column)
    @test gf2_record.metadata.source == :supplied_link_witness
    @test gf2_record.verification.overall_ok == true
    @test gf2_record.verification.metadata_ok == true
    @test gf2_record.verification.selected_monic_ok == true
    @test gf2_record.verification.tail_reduction_ok == true
    @test gf2_record.verification.resultants_ok == true
    @test gf2_record.verification.bezout_ok == true
    @test gf2_record.verification.coverage_ok == true
    @test gf2_record.verification.path_ok == true
    @test Suslin.verify_ecp_link_witness(gf2_record) == true

    qq_entry = _case_by_id("ecp-monic-first-entry-qq")
    qq_column = _column(qq_entry)
    qq_witness = _qq_link_witness(qq_entry)
    qq_record = Suslin.ecp_link_witness(
        qq_column,
        qq_entry.ring.object;
        variable_order = qq_entry.ring.generators,
        selected_variable = qq_entry.ring.generators[1],
        supplied_link_witness = qq_witness,
    )
    @test !any(is_unit, qq_column)
    @test qq_record.metadata.source == :supplied_link_witness
    @test qq_record.verification.overall_ok == true
    @test qq_record.verification.metadata_ok == true
    @test qq_record.verification.selected_monic_ok == true
    @test qq_record.verification.lengths_ok == true
    @test qq_record.verification.tail_reduction_ok == true
    @test qq_record.verification.resultants_ok == true
    @test qq_record.verification.bezout_ok == true
    @test qq_record.verification.coverage_ok == true
    @test qq_record.verification.path_ok == true
    @test Suslin.verify_ecp_link_witness(qq_record) == true

    x, y = qq_entry.ring.generators
    R = qq_entry.ring.object
    @test_throws ArgumentError Suslin.ecp_link_witness(
        qq_column,
        R;
        variable_order = qq_entry.ring.generators,
        selected_variable = x,
        supplied_link_witness = _mutate_witness(
            qq_witness;
            resultants = (y^2 + one(R), y + one(R)),
        ),
    )
    @test_throws ArgumentError Suslin.ecp_link_witness(
        qq_column,
        R;
        variable_order = qq_entry.ring.generators,
        selected_variable = x,
        supplied_link_witness = _mutate_witness(
            qq_witness;
            bezout_coefficients = ((; f = one(R), h = y), qq_witness.bezout_coefficients[2]),
        ),
    )
    @test_throws ArgumentError Suslin.ecp_link_witness(
        qq_column,
        R;
        variable_order = qq_entry.ring.generators,
        selected_variable = x,
        supplied_link_witness = _mutate_witness(
            qq_witness;
            coverage_multipliers = (one(R), one(R)),
        ),
    )
    @test_throws ArgumentError Suslin.ecp_link_witness(
        qq_column,
        R;
        variable_order = qq_entry.ring.generators,
        selected_variable = x,
        supplied_link_witness = _mutate_witness(
            qq_witness;
            path_points = (zero(R), y^2 * x + one(R), x),
        ),
    )

    bad_verification = merge(qq_record.verification, (; coverage_ok = false, overall_ok = false))
    @test !Suslin.verify_ecp_link_witness(_tamper_record_field(qq_record, :verification, bad_verification))
    @test !Suslin.verify_ecp_link_witness(_tamper_record_field(qq_record, :path_points, _replace_tuple_entry(qq_record.path_points, 2, y^2 * x + one(R))))
end

@testset "ECP link witnesses remain staged without supplied metadata" begin
    entry = _case_by_id("ecp-monic-first-entry-qq")
    column = _column(entry)
    R = entry.ring.object
    x = entry.ring.generators[1]
    err = try
        Suslin.ecp_link_witness(
            column,
            R;
            variable_order = entry.ring.generators,
            selected_variable = x,
        )
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("supplied_link_witness", sprint(showerror, err))
end
