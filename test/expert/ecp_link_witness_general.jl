using Test
using Oscar
using Suslin

include(joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl"))

function _lw_column(entry)
    return [getproperty(entry.column_entries, name) for name in entry.column_order]
end

function _lw_context(entry)
    return Suslin.ecp_input_context(
        _lw_column(entry),
        entry.ring.object;
        variable_order = entry.ring.generators,
        selected_variable = entry.selected_variable.generator,
        unimodularity_witness = entry.unimodularity.coefficients,
    )
end

function _lw_normalization(entry)
    return Suslin.ecp_monicity_normalization(
        _lw_context(entry);
        selected_variable = entry.selected_variable.generator,
        max_shift_power = 2,
    )
end

function _lw_recompute_tail(record, idx::Int)
    total = zero(record.ring)
    tail_entries = record.original_column[2:end]
    for tail_idx in eachindex(tail_entries)
        total += record.tail_reductions[idx].lifted_tail_coefficients[tail_idx] * tail_entries[tail_idx]
    end
    return total
end

function _lw_replace_tuple_entry(values::Tuple, idx::Int, value)
    return ntuple(j -> j == idx ? value : values[j], length(values))
end

function _lw_replace_record_field(record, field::Symbol, value)
    fields = fieldnames(typeof(record))
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown ECP link witness field $(field)")
    values = [getfield(record, name) for name in fields]
    values[idx] = value
    return typeof(record)(values...)
end

function _lw_assert_replayed_equations(record)
    @test record isa Suslin.ECPLinkWitnessRecord
    @test record.metadata.source == :extracted_link_witness
    @test Suslin.verify_ecp_link_witness(record)
    @test record.verification.tail_reduction_ok
    @test record.verification.resultants_ok
    @test record.verification.bezout_ok
    @test record.verification.coverage_ok
    @test record.verification.path_ok

    for idx in eachindex(record.tail_reductions)
        G = _lw_recompute_tail(record, idx)
        @test record.tail_reductions[idx].G == G
        @test record.tail_reductions[idx].tilde_G == G
        @test record.resultants[idx] == resultant(record.selected_monic_entry, G, record.selected_variable_index)
        bezout = record.bezout_coefficients[idx]
        @test bezout.f * record.selected_monic_entry + bezout.h * G == record.resultants[idx]
    end

    coverage_total = sum(
        record.coverage_multipliers[idx] * record.resultants[idx] for idx in eachindex(record.resultants);
        init = zero(record.ring),
    )
    @test coverage_total == one(record.ring)
    @test first(record.path_points) == zero(record.ring)
    @test last(record.path_points) == record.selected_variable
    for idx in eachindex(record.resultants)
        @test record.path_points[idx + 1] - record.path_points[idx] ==
              record.coverage_multipliers[idx] * record.resultants[idx] * record.selected_variable
    end
end

@testset "ECP link witness extraction replay and diagnostics" begin
    cases = ECPMainlineFixtureCatalog.cases_by_id()

    qq_entry = cases["ecp-mainline-qq-link-bezout"]
    qq_record = Suslin.ecp_link_witness(_lw_normalization(qq_entry))
    _lw_assert_replayed_equations(qq_record)
    @test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(qq_record, :resultants, _lw_replace_tuple_entry(qq_record.resultants, 1, qq_record.resultants[1] + one(qq_record.ring))))
    @test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(qq_record, :bezout_coefficients, _lw_replace_tuple_entry(qq_record.bezout_coefficients, 1, merge(qq_record.bezout_coefficients[1], (; f = qq_record.bezout_coefficients[1].f + one(qq_record.ring))))))
    @test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(qq_record, :path_points, _lw_replace_tuple_entry(qq_record.path_points, 2, qq_record.path_points[2] + one(qq_record.ring))))
    @test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(qq_record, :tail_reductions, _lw_replace_tuple_entry(qq_record.tail_reductions, 1, merge(qq_record.tail_reductions[1], (; lifted_tail_coefficients = _lw_replace_tuple_entry(qq_record.tail_reductions[1].lifted_tail_coefficients, 1, qq_record.tail_reductions[1].lifted_tail_coefficients[1] + one(qq_record.ring)))))))

    length4_entry = cases["ecp-mainline-length4-coupled-qq"]
    length4_record = Suslin.ecp_link_witness(_lw_normalization(length4_entry))
    _lw_assert_replayed_equations(length4_record)
    failure = Suslin.ecp_link_witness(_lw_normalization(length4_entry); max_tail_terms = 1)
    @test failure isa Suslin.ECPLinkWitnessExtractionFailure
    @test failure.kind == :link_witness_cover_not_proved
end
