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
    lazy_metadata = Suslin._factor_laurent_gl_lazy_determinant_peel(
        A;
        progress_callback = record -> push!(lazy_progress, record),
        determinant_probe = _issue155_lazy_probe(original_size[1], lazy_progress, lazy_probes),
    )

    @test lazy_metadata.determinant_source == :deferred_submatrix
    @test lazy_metadata.determinant_classification == :laurent_monomial_unit
    @test lazy_metadata.supported
    @test lazy_metadata.normalized_deferred_core !== nothing
    @test det(lazy_metadata.normalized_deferred_core) ==
        one(base_ring(lazy_metadata.normalized_deferred_core))
    @test lazy_metadata.staged_boundary === nothing
    @test !isempty(lazy_progress)
    @test any(record -> record.completed_steps >= 1, lazy_progress)
    first_completed_progress = first(record for record in lazy_progress if record.completed_steps >= 1)
    @test first_completed_progress.current_dimension < original_size[1]
    @test first_completed_progress.last_completed_dimension == original_size[1]
    @test !isempty(lazy_probes)
    @test first(lazy_probes).size[1] < original_size[1]
    @test first(lazy_probes).size[2] < original_size[2]
    @test first(lazy_probes).completed_before_probe >= 1

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

@testset "lazy Laurent peel determinant-one continuation" begin
    entry = _issue155_fixture("determinant-one-triangular")
    A = entry.inputs.matrix
    original_size = (nrows(A), ncols(A))

    progress_records = Any[]
    probe_records = Any[]
    certificate = Suslin._factor_laurent_gl_lazy_determinant_peel(
        A;
        progress_callback = record -> push!(progress_records, record),
        determinant_probe = candidate -> begin
            profile = Suslin.classify_laurent_determinant(candidate)
            push!(probe_records, (;
                size = (nrows(candidate), ncols(candidate)),
                classification = profile.classification,
                completed_before_probe = _issue155_max_completed_steps(progress_records),
            ))
            return profile
        end,
    )

    @test Suslin._verify_laurent_column_peel_replay(certificate)
    @test verify_factorization(A, certificate.factors)
    @test length(certificate.peel_steps) == original_size[1] - 2
    @test !isempty(probe_records)
    @test first(probe_records).classification == :one
    @test first(probe_records).size[1] < original_size[1]
    @test first(probe_records).size[2] < original_size[2]
    @test first(probe_records).completed_before_probe >= 1
end

@testset "lazy Laurent peel size guard" begin
    entry = _issue155_fixture("determinant-one-triangular")
    R = entry.ring.object
    too_small = identity_matrix(R, 2)

    err = try
        Suslin._factor_laurent_gl_lazy_determinant_peel(too_small)
        nothing
    catch caught
        caught
    end

    @test err isa ArgumentError
    @test occursin("requires size at least 3", sprint(showerror, err))
end
