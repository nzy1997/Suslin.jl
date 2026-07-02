using Test

const POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _normality_boundary_read(path_parts...)
    return read(joinpath(POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT, path_parts...), String)
end

function _normality_boundary_required_phrases()
    return [
        "ordinary-polynomial normality/conjugation certificates",
        "Cohn-type realization certificates",
        "rank-one normality certificates",
        "conjugated-elementary normality certificates",
        "staged ECP induction/normality adapter replays a nested conjugated-elementary certificate",
        "Murthy local `SL_3` solver (#182) is supported for the proven ordinary/local-witness contract",
        "ordinary factor vectors are exposed only when the certificate can materialize them over the base ring",
        "ordinary-polynomial ECP unimodular-column reducer (#185) is accepted",
        "Polynomial column peel records the verified ECP certificate used for each last-column peel step",
        "Quillen automatic patching (#183), general `SL_3` (#184), recursive `SL_n` (#186)",
        "full public Park-Woodburn acceptance (#187), coefficient-ring support beyond exact field-backed ordinary polynomial rings",
        "Laurent/ToricBuilder mainline acceptance, and Steinberg factor-count optimization remain staged boundaries",
    ]
end

function _normality_boundary_contains(text::AbstractString, phrase::AbstractString)
    normalized_text = replace(text, r"\s+" => " ")
    normalized_phrase = replace(phrase, r"\s+" => " ")
    return occursin(normalized_phrase, normalized_text)
end

@testset "polynomial normality support boundary documentation" begin
    docs = Dict(
        "README.md" => _normality_boundary_read("README.md"),
        "docs/src/index.md" => _normality_boundary_read("docs", "src", "index.md"),
    )

    for (path, text) in docs
        @testset "$path" begin
            for phrase in _normality_boundary_required_phrases()
                @test _normality_boundary_contains(text, phrase)
            end
            @test !_normality_boundary_contains(text, "full #187 acceptance is supported")
            @test !_normality_boundary_contains(text, "ToricBuilder `case_008` is supported")
            @test !_normality_boundary_contains(text, "Steinberg optimization is supported")
            @test !_normality_boundary_contains(text, "Laurent `GL_n` is fully supported")
        end
    end
end

@testset "polynomial normality certificate expert gate registration" begin
    runtests = _normality_boundary_read("test", "runtests.jl")
    for expert_file in (
        "expert/cohn_type.jl",
        "expert/normality_rank_one.jl",
        "expert/normality.jl",
        "expert/ecp_induction_normality.jl",
    )
        @test occursin("\"$(expert_file)\"", runtests)
    end
    @test occursin("\"expert/polynomial_normality_support_boundary.jl\"", runtests)
end
