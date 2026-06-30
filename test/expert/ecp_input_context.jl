using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl"))

struct _InputContextZeroBasedVector{T} <: AbstractVector{T}
    data::Vector{T}
end

Base.size(v::_InputContextZeroBasedVector) = (length(v.data),)
Base.axes(v::_InputContextZeroBasedVector) = (0:(length(v.data) - 1),)
Base.firstindex(v::_InputContextZeroBasedVector) = 0
Base.lastindex(v::_InputContextZeroBasedVector) = length(v.data) - 1
Base.getindex(v::_InputContextZeroBasedVector, idx::Int) = v.data[idx + 1]

struct _InputContextThrowingLengthVector{T} <: AbstractVector{T}
    data::Vector{T}
end

Base.size(v::_InputContextThrowingLengthVector) = (length(v.data),)
Base.axes(v::_InputContextThrowingLengthVector) = (Base.OneTo(length(v.data)),)
Base.firstindex(::_InputContextThrowingLengthVector) = 1
Base.lastindex(v::_InputContextThrowingLengthVector) = length(v.data)
Base.length(::_InputContextThrowingLengthVector) = throw(ErrorException("length unavailable"))
Base.getindex(v::_InputContextThrowingLengthVector, idx::Int) = v.data[idx]

struct _InputContextThrowingGetindexVector{T} <: AbstractVector{T}
    data::Vector{T}
end

Base.size(v::_InputContextThrowingGetindexVector) = (length(v.data),)
Base.axes(v::_InputContextThrowingGetindexVector) = (Base.OneTo(length(v.data)),)
Base.getindex(::_InputContextThrowingGetindexVector, ::Int) =
    throw(ErrorException("entry unavailable"))

struct _InputContextBadInteger <: Integer end

Base.Int(::_InputContextBadInteger) = throw(ErrorException("index unavailable"))

mutable struct _InputContextSecondPassOrder
    data::Tuple
    passes::Int
end

function Base.iterate(order::_InputContextSecondPassOrder, state::Int = 1)
    if state == 1
        order.passes += 1
        order.passes > 1 && throw(ErrorException("variable order unavailable"))
    end
    state > length(order.data) && return nothing
    return order.data[state], state + 1
end

Base.:(==)(order::_InputContextSecondPassOrder, rhs::Tuple) = order.data == rhs

struct _InputContextRingKindThrowRing end

Suslin._is_laurent_polynomial_ring(::_InputContextRingKindThrowRing) =
    throw(ErrorException("ring kind unavailable"))

struct _InputContextWitnessThrowRing end

AbstractAlgebra.gens(::_InputContextWitnessThrowRing) = (:x,)
Base.zero(::_InputContextWitnessThrowRing) = 0
Base.one(::_InputContextWitnessThrowRing) = 1
(::_InputContextWitnessThrowRing)(value) = value
Suslin._is_laurent_polynomial_ring(::_InputContextWitnessThrowRing) = false
Suslin.is_unimodular_column(::AbstractVector, ::_InputContextWitnessThrowRing) = true
Suslin._unimodular_witness(::AbstractVector, ::_InputContextWitnessThrowRing) =
    throw(ErrorException("witness unavailable"))

struct _InputContextOneThrowRing end

AbstractAlgebra.gens(::_InputContextOneThrowRing) = (:x,)
Base.zero(::_InputContextOneThrowRing) = 0
Base.one(::_InputContextOneThrowRing) = throw(ErrorException("unit unavailable"))
(::_InputContextOneThrowRing)(value) = value
Suslin._is_laurent_polynomial_ring(::_InputContextOneThrowRing) = false

function _input_context_fixture_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _input_context_mainline_column(entry)
    return [getproperty(entry.column_entries, name) for name in entry.column_order]
end

function _input_context_witness_total(ctx)
    total = zero(ctx.ring)
    for idx in 1:ctx.column_length
        total += ctx.unimodularity_witness[idx] * ctx.column[idx]
    end
    return total
end

function _input_context_staged_failure_reason(diagnostic)
    return diagnostic.status == :unsupported ? diagnostic.failure_code : nothing
end

function _input_context_replace_field(ctx, field::Symbol, value)
    fields = fieldnames(typeof(ctx))
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown ECPInputContext field $(field)")
    values = [getfield(ctx, name) for name in fields]
    values[idx] = value
    return Suslin.ECPInputContext(values...)
end

function _assert_checked_context(
    ctx,
    column,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
)
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    normalized_order = tuple(Suslin._ecp_normalize_variable_order(R, variable_order)...)

    @test ctx isa Suslin.ECPInputContext
    @test Suslin.verify_ecp_input_context(ctx)
    @test ctx.column == [Suslin._coerce_into_ring(R, column[idx], "column[$idx]") for idx in eachindex(column)]
    @test ctx.ring == R
    @test ctx.ring_profile == Suslin._column_reduction_ring_profile(R)
    @test ctx.variables == tuple(gens(R)...)
    @test ctx.variable_order == normalized_order
    @test ctx.column_length == length(column)
    @test _input_context_witness_total(ctx) == one(R)
    @test ctx.unimodularity_witness == Suslin._unimodular_witness(ctx.column, R)
    if selected_variable === nothing
        @test ctx.selected_variable_index === nothing
        @test ctx.selected_variable === nothing
    else
        selected_index = Suslin._ecp_selected_variable_index(R, selected_variable)
        @test ctx.selected_variable_index == selected_index
        @test ctx.selected_variable == gens(R)[selected_index]
        @test count(==(ctx.selected_variable), ctx.variable_order) == 1
    end
    @test ctx.support_classification == diagnostic.status
    @test ctx.staged_failure_reason == _input_context_staged_failure_reason(diagnostic)
    @test ctx.staged_diagnostic == diagnostic
    @test ctx.verification.overall_ok
    @test ctx.verification.one_based_indexing_ok
    @test ctx.verification.column_ok
    @test ctx.verification.ring_profile_ok
    @test ctx.verification.variables_ok
    @test ctx.verification.variable_order_ok
    @test ctx.verification.column_length_ok
    @test ctx.verification.unimodularity_witness_ok
    @test ctx.verification.selected_variable_ok
    @test ctx.verification.staged_diagnostic_ok
end

@testset "checked ECP input contexts" begin
    mainline_cases = ECPMainlineFixtureCatalog.cases_by_id()
    column_cases = ECPColumnFixtureCatalog.cases_by_id()

    supported_entry = mainline_cases["ecp-mainline-gf2-hard-slice"]
    supported_R = supported_entry.ring.object
    supported_column = _input_context_mainline_column(supported_entry)
    supported_ctx = Suslin.ecp_input_context(
        supported_column,
        supported_R;
        variable_order = supported_entry.ring.generators,
        selected_variable = supported_entry.selected_variable.generator,
        unimodularity_witness = supported_entry.unimodularity.coefficients,
    )
    _assert_checked_context(
        supported_ctx,
        supported_column,
        supported_R;
        variable_order = supported_entry.ring.generators,
        selected_variable = supported_entry.selected_variable.generator,
    )
    @test supported_ctx.support_classification == :supported
    @test supported_ctx.staged_failure_reason === nothing

    staged_entry = column_cases["ecp-unsupported-unimodular-gf2"]
    staged_R = staged_entry.ring.object
    staged_column = _input_context_fixture_column(staged_entry)
    staged_ctx = Suslin.ecp_input_context(
        staged_column,
        staged_R;
        variable_order = tuple(gens(staged_R)...),
        selected_variable = gens(staged_R)[1],
    )
    _assert_checked_context(
        staged_ctx,
        staged_column,
        staged_R;
        selected_variable = gens(staged_R)[1],
    )
    @test staged_ctx.support_classification == :unsupported
    @test staged_ctx.staged_failure_reason == :unsupported_polynomial_column_family

    default_selected_ctx = Suslin.ecp_input_context(staged_column, staged_R)
    _assert_checked_context(default_selected_ctx, staged_column, staged_R)
    @test default_selected_ctx.selected_variable_index === nothing
    @test default_selected_ctx.selected_variable === nothing

    unsupported_column = [
        zero(staged_R),
        gens(staged_R)[1]^2,
        gens(staged_R)[1] * gens(staged_R)[2] + one(staged_R),
        zero(staged_R),
    ]
    @test Suslin.is_unimodular_column(unsupported_column, staged_R)
    unsupported_ctx = Suslin.ecp_input_context(
        unsupported_column,
        staged_R;
        selected_variable = gens(staged_R)[1],
    )
    _assert_checked_context(
        unsupported_ctx,
        unsupported_column,
        staged_R;
        selected_variable = gens(staged_R)[1],
    )
    @test unsupported_ctx.column_length == 4
    @test unsupported_ctx.support_classification == :unsupported
    @test unsupported_ctx.staged_failure_reason == :unsupported_polynomial_column_family

    length4_entry = mainline_cases["ecp-mainline-length4-coupled-qq"]
    length4_R = length4_entry.ring.object
    length4_column = _input_context_mainline_column(length4_entry)
    length4_ctx = Suslin.ECPInputContext(
        length4_column,
        length4_R;
        variable_order = length4_entry.ring.generators,
        selected_variable = length4_entry.selected_variable.generator,
        unimodularity_witness = length4_entry.unimodularity.coefficients,
    )
    _assert_checked_context(
        length4_ctx,
        length4_column,
        length4_R;
        variable_order = length4_entry.ring.generators,
        selected_variable = length4_entry.selected_variable.generator,
    )
    length4_diagnostic = Suslin.diagnose_unimodular_column_reduction(length4_column, length4_R)
    @test length4_ctx.support_classification == length4_diagnostic.status
    @test length4_ctx.staged_failure_reason == _input_context_staged_failure_reason(length4_diagnostic)

    tampered_witness = collect(length4_ctx.unimodularity_witness)
    tampered_witness[1] += one(length4_R)
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :unimodularity_witness, tampered_witness),
    )

    @test !Suslin.verify_ecp_input_context(nothing)
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :column,
            _InputContextZeroBasedVector(collect(length4_ctx.column)),
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :column,
            _InputContextThrowingLengthVector(collect(length4_ctx.column)),
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :column,
            _InputContextThrowingGetindexVector(collect(length4_ctx.column)),
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :column, typeof(length4_ctx.column)()),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :ring,
            _InputContextRingKindThrowRing(),
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :ring, _InputContextWitnessThrowRing()),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :ring, _InputContextOneThrowRing()),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :ring, nothing),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :ring_profile, (; tampered = true)),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :variables, reverse(length4_ctx.variables)),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :variable_order, (:not_a_generator,)),
    )
    throwing_order = _InputContextSecondPassOrder(length4_ctx.variable_order, 0)
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :variable_order, throwing_order),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :column_length, length4_ctx.column_length + 1),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :unimodularity_witness,
            _InputContextZeroBasedVector(collect(length4_ctx.unimodularity_witness)),
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :unimodularity_witness,
            _InputContextThrowingLengthVector(collect(length4_ctx.unimodularity_witness)),
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :unimodularity_witness,
            length4_ctx.unimodularity_witness[1:(end - 1)],
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :selected_variable_index,
            _InputContextBadInteger(),
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :selected_variable,
            length4_ctx.variables[1] + length4_ctx.variables[2],
        ),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :support_classification, :tampered),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :staged_failure_reason, :tampered),
    )
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(
            length4_ctx,
            :staged_diagnostic,
            merge(length4_ctx.staged_diagnostic, (; message = "tampered")),
        ),
    )
end

@testset "ECP input context rejects invalid boundaries" begin
    column_cases = ECPColumnFixtureCatalog.cases_by_id()
    non_unimodular_entry = column_cases["ecp-non-unimodular-gf2"]
    non_unimodular_R = non_unimodular_entry.ring.object
    non_unimodular = _input_context_fixture_column(non_unimodular_entry)
    @test_throws ArgumentError Suslin.ecp_input_context(non_unimodular, non_unimodular_R)

    staged_entry = column_cases["ecp-unsupported-unimodular-gf2"]
    R = staged_entry.ring.object
    x, y = gens(R)
    staged_column = _input_context_fixture_column(staged_entry)
    @test_throws ArgumentError Suslin.ecp_input_context([one(R), x], R)
    @test_throws ArgumentError Suslin.ecp_input_context(staged_column, R; selected_variable = x + y)
    @test_throws ArgumentError Suslin.ecp_input_context(
        staged_column,
        R;
        variable_order = (x,),
        selected_variable = y,
    )

    mainline_cases = ECPMainlineFixtureCatalog.cases_by_id()
    length4_entry = mainline_cases["ecp-mainline-length4-coupled-qq"]
    length4_R = length4_entry.ring.object
    length4_column = _input_context_mainline_column(length4_entry)
    bad_witness = collect(length4_entry.unimodularity.coefficients)
    bad_witness[1] = zero(length4_R)
    @test_throws ArgumentError Suslin.ecp_input_context(
        length4_column,
        length4_R;
        selected_variable = length4_entry.selected_variable.generator,
        unimodularity_witness = bad_witness,
    )
    @test_throws ArgumentError Suslin.ecp_input_context(
        length4_column,
        length4_R;
        selected_variable = length4_entry.selected_variable.generator,
        unimodularity_witness = _InputContextZeroBasedVector(
            collect(length4_entry.unimodularity.coefficients),
        ),
    )

    laurent_R, (lx, ly) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    laurent_column = [
        lx^-1 * ly^-1 * (lx + ly^2),
        lx^-1 * ly^-1 * (lx * ly + lx + one(laurent_R)),
        lx^-1 * ly^-1 * (lx^2 + lx * ly + ly + one(laurent_R)),
    ]
    @test_throws ArgumentError Suslin.ecp_input_context(laurent_column, laurent_R)
end
