using Test
using Suslin
using Oscar

const SUPPORT_BOUNDARY_EVIDENCE_PAGE = joinpath(
    @__DIR__,
    "..",
    "..",
    "docs",
    "audits",
    "2026-06-26-laurent-toricbuilder-support-boundary-evidence.md",
)

function _read_support_boundary_evidence()
    @test isfile(SUPPORT_BOUNDARY_EVIDENCE_PAGE)
    return read(SUPPORT_BOUNDARY_EVIDENCE_PAGE, String)
end

@testset "documentation smoke" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = matrix(R, [
        one(R)      one(R) + X      zero(R);
        X           one(R) + X + X^2 zero(R);
        zero(R)     zero(R)         one(R)
    ])

    factors = elementary_factorization(A)
    @test verify_factorization(A, factors)

    @testset "support boundary evidence page" begin
        evidence = _read_support_boundary_evidence()
        @test occursin(
            "| `case_010` ToricBuilder Q-block | `laurent_gl_factorization_certificate` | `gl_certificate_pass`; verified `true`; decomposed base matrices `48` | public `elementary_factorization` remains `staged_boundary` for the original Laurent `GL_n` input |",
            evidence,
        )
        @test occursin(
            "| `case_008` bounded exercise | bounded Laurent `GL_n` certificate route | `certified_algorithm_boundary` at `certificate_construction` under explicit `--exercise=case_008 --timeout-seconds=120` | not a default report pass; remains a staged algorithm boundary |",
            evidence,
        )
        @test occursin(
            "ToricBuilder Q-block -> classify Laurent determinant -> normalize Laurent GL_n determinant -> factor determinant-one core -> verify Laurent GL_n certificate",
            evidence,
        )
        @test occursin(
            "elementary_factorization(A) -> exact elementary factor sequence -> verify_factorization(A, factors)",
            evidence,
        )
        @test occursin(
            "original-input `elementary_factorization` for Laurent `GL_n` remains a `staged boundary`",
            evidence,
        )
        @test occursin("not arbitrary Park-Woodburn", evidence)
        @test occursin("not arbitrary Laurent `GL_n`", evidence)
        @test occursin("julia --project=. test/runtests.jl", evidence)
        @test occursin("julia --project=. test/runtests.jl expert", evidence)
    end
end
