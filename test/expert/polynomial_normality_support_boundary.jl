using Test

const POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _normality_boundary_read(path_parts...)
    return read(joinpath(POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT, path_parts...), String)
end

function _normality_boundary_required_phrases()
    return [
        "ordinary-polynomial normality/conjugation certificates",
        "ordinary-polynomial ECP unimodular-column reducer (#185) is accepted",
        "recursive ordinary-polynomial `SL_n` driver (#186) is supported for exact field-backed ordinary-polynomial `SL_n`, `n > 3`, inputs",
        "whose recursive peel steps verify #185 ECP evidence",
        "whose final `SL_3` block verifies #184 route evidence",
        "public route certificate carries #186 mainline provenance",
        "determinant-one `SL_n` inputs missing ECP, final `SL_3`, variable, local-form, Quillen/Murthy, or unsupported-ring evidence remain staged",
        "`:missing_ecp_evidence`",
        "`:missing_final_sl3_route`",
        "legacy fast-local/disjoint-block examples may still verify factors but do not count as #186 mainline support by themselves",
        "full public Park-Woodburn acceptance (#187), coefficient-ring support beyond exact field-backed ordinary polynomial rings",
        "Laurent/ToricBuilder mainline acceptance, and Steinberg factor-count optimization remain staged boundaries",
    ]
end

function _normality_boundary_contains(text::AbstractString, phrase::AbstractString)
    normalized_text = replace(text, r"\s+" => " ")
    normalized_phrase = replace(phrase, r"\s+" => " ")
    return occursin(normalized_phrase, normalized_text)
end

function _normality_boundary_audit_path()
    return joinpath(
        POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT,
        "docs",
        "audits",
        "2026-07-03-issue-186-recursive-sln-acceptance.md",
    )
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
            @test !_normality_boundary_contains(text, "#187 final mainline acceptance is supported")
            @test !_normality_boundary_contains(text, "arbitrary Laurent `GL_n` support is complete")
            @test !_normality_boundary_contains(text, "ToricBuilder support is complete")
            @test !_normality_boundary_contains(text, "Steinberg factor-count optimization is supported")
        end
    end
end

@testset "issue 186 recursive SLn acceptance audit note" begin
    audit_path = _normality_boundary_audit_path()
    @test isfile(audit_path)
    audit = read(audit_path, String)
    for phrase in (
            "#260",
            "#261",
            "#262",
            "#263",
            "#264",
            "#265",
            "#266",
            "Park-Woodburn Section 3",
            "Park-Woodburn Section 4",
            "Park-Woodburn Section 5",
            ":missing_ecp_evidence",
            ":missing_final_sl3_route",
            "does not close #187")
        @test _normality_boundary_contains(audit, phrase)
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
