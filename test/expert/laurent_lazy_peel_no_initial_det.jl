using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

struct _Issue155InitialDeterminantProbeError <: Exception
    message::String
end

Base.showerror(io::IO, err::_Issue155InitialDeterminantProbeError) =
    print(io, err.message)

function _issue155_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue155_max_completed_steps(progress_records)
    isempty(progress_records) && return 0
    return maximum(record.completed_steps for record in progress_records)
end

function _issue155_lazy_probe(original_dimension::Int, progress_records, probe_records)
    return function (candidate)
        completed_before_probe = _issue155_max_completed_steps(progress_records)
        push!(probe_records, (;
            size = (nrows(candidate), ncols(candidate)),
            completed_before_probe,
        ))
        if nrows(candidate) == original_dimension || completed_before_probe < 1
            throw(_Issue155InitialDeterminantProbeError(
                "initial determinant classification invoked before lazy Laurent peel completed a step",
            ))
        end
        return Suslin.classify_laurent_determinant(candidate)
    end
end

function _issue155_eager_probe(probe_records)
    return function (candidate)
        push!(probe_records, (; size = (nrows(candidate), ncols(candidate))))
        throw(_Issue155InitialDeterminantProbeError(
            "eager determinant classification invoked at the original Laurent matrix",
        ))
    end
end

@testset "lazy Laurent peel defers initial determinant classification" begin
    entry = _issue155_fixture("monomial-unit-row-column-cores")
    A = entry.inputs.matrix
    original_size = (nrows(A), ncols(A))

    lazy_progress = Any[]
    lazy_probes = Any[]
    lazy_err = try
        Suslin._factor_laurent_gl_lazy_determinant_peel(
            A;
            progress_callback = record -> push!(lazy_progress, record),
            determinant_probe = _issue155_lazy_probe(original_size[1], lazy_progress, lazy_probes),
        )
        nothing
    catch err
        err
    end

    @test lazy_err isa ArgumentError
    @test occursin("lazy Laurent determinant correction", sprint(showerror, lazy_err))
    @test !isempty(lazy_progress)
    @test any(record -> record.completed_steps >= 1, lazy_progress)
    first_completed_progress = first(record for record in lazy_progress if record.completed_steps >= 1)
    @test first_completed_progress.current_dimension < original_size[1]
    @test first_completed_progress.last_completed_dimension == original_size[1]
    @test !isempty(lazy_probes)
    @test first(lazy_probes).size[1] < original_size[1]
    @test first(lazy_probes).size[2] < original_size[2]
    @test first(lazy_probes).completed_before_probe >= 1
    @test !(lazy_err isa _Issue155InitialDeterminantProbeError)

    eager_progress = Any[]
    eager_probes = Any[]
    eager_err = try
        Suslin._factor_laurent_sl_column_peel(
            A;
            progress_callback = record -> push!(eager_progress, record),
            determinant_probe = _issue155_eager_probe(eager_probes),
        )
        nothing
    catch err
        err
    end

    @test eager_err isa _Issue155InitialDeterminantProbeError
    @test occursin("eager determinant classification", sprint(showerror, eager_err))
    @test !isempty(eager_progress)
    @test _issue155_max_completed_steps(eager_progress) == 0
    @test !isempty(eager_probes)
    @test first(eager_probes).size == original_size
end
