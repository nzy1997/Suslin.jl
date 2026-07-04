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
const ISSUE187_ACCEPTANCE_AUDIT_PAGE = joinpath(
    @__DIR__,
    "..",
    "..",
    "docs",
    "audits",
    "2026-07-04-issue-187-park-woodburn-mainline-acceptance.md",
)
const ISSUE188_ACCEPTANCE_AUDIT_PAGE = joinpath(
    @__DIR__,
    "..",
    "..",
    "docs",
    "audits",
    "2026-07-04-issue-188-steinberg-optimization.md",
)
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const README_PATH = joinpath(REPO_ROOT, "README.md")
const DOCS_INDEX_PATH = joinpath(REPO_ROOT, "docs", "src", "index.md")

function _read_support_boundary_evidence()
    @test isfile(SUPPORT_BOUNDARY_EVIDENCE_PAGE)
    return read(SUPPORT_BOUNDARY_EVIDENCE_PAGE, String)
end

function _read_issue187_acceptance_audit()
    @test isfile(ISSUE187_ACCEPTANCE_AUDIT_PAGE)
    return read(ISSUE187_ACCEPTANCE_AUDIT_PAGE, String)
end

function _read_issue188_acceptance_audit()
    @test isfile(ISSUE188_ACCEPTANCE_AUDIT_PAGE)
    return read(ISSUE188_ACCEPTANCE_AUDIT_PAGE, String)
end

function _read_repo_text(path)
    @test isfile(path)
    return read(path, String)
end

function _paragraphs(text)
    return split(replace(text, "\r\n" => "\n"), "\n\n")
end

function _squash_whitespace(text)
    return replace(strip(text), r"\s+" => " ")
end

function _assert_not_claimed_as_issue187(text, item)
    for paragraph in _paragraphs(text)
        squashed = _squash_whitespace(paragraph)
        lower_squashed = lowercase(squashed)
        if occursin("#187", lower_squashed) && occursin(lowercase(item), lower_squashed)
            @test occursin("separate", squashed) ||
                  occursin("out of scope", squashed) ||
                  occursin("not part of", squashed) ||
                  occursin("outside", lowercase(squashed))
        end
    end
end

function _assert_issue187_public_contract(text)
    squashed = _squash_whitespace(text)
    @test occursin(
        "The final ordinary-polynomial Park-Woodburn public contract (#187) is supported",
        squashed,
    )
    @test occursin(
        "exact field-backed ordinary-polynomial determinant-one `SL_3`",
        squashed,
    )
    @test occursin(
        "exact field-backed ordinary-polynomial determinant-one `SL_n`, `n > 3`",
        squashed,
    )
    @test occursin("implemented evidence-backed route", squashed)
    @test occursin(
        "Unsupported coefficient rings remain out of scope",
        squashed,
    )
    @test occursin(
        "Arbitrary Laurent `GL_n`, ToricBuilder mainline acceptance, and Steinberg factor-count optimization (#188) remain separate from #187",
        squashed,
    )
    @test !occursin("full public Park-Woodburn acceptance (#187)", squashed)
    _assert_not_claimed_as_issue187(text, "unsupported coefficient rings")
    _assert_not_claimed_as_issue187(text, "arbitrary Laurent `GL_n`")
    _assert_not_claimed_as_issue187(text, "ToricBuilder")
    _assert_not_claimed_as_issue187(text, "factor-count optimization")
    @test occursin("#187 closeout coverage audit", squashed)
end

function _assert_issue187_acceptance_audit(text)
    squashed = _squash_whitespace(text)
    for pair in (
        ("#181", "#195"),
        ("#182", "#212"),
        ("#183", "#220"),
        ("#184", "#239"),
        ("#185", "#249"),
        ("#186", "#266"),
    )
        @test occursin(pair[1], squashed)
        @test occursin(pair[2], squashed)
    end
    for issue_id in ("#270", "#271", "#272")
        @test occursin(issue_id, squashed)
    end
    for case_id in (
        "pw-mainline-sl3-multivariate-issue184-qq",
        "pw-mainline-sln-recursive-issue185-186-gf2",
        "pw-mainline-readme-ordinary-polynomial-qq",
        "pw-mainline-negative-unsupported-coefficient-ring",
        "pw-mainline-negative-missing-ecp-evidence",
        "pw-mainline-negative-missing-final-sl3-evidence",
        "pw-mainline-negative-laurent-boundary",
    )
        @test occursin(case_id, squashed)
    end
    @test occursin(
        "Steinberg factor-count optimization (#188) remains separate",
        squashed,
    )
    @test occursin(
        "Laurent/ToricBuilder mainline support remains separate",
        squashed,
    )
    @test occursin("Unsupported coefficient rings remain negative controls", squashed)
end

function _assert_issue188_no_overclaims(text)
    for paragraph in _paragraphs(text)
        squashed = _squash_whitespace(paragraph)
        lower_squashed = lowercase(squashed)
        mentions_issue188 =
            occursin("#188", lower_squashed) || occursin("steinberg", lower_squashed)

        if mentions_issue188 && occursin("enabled by default", lower_squashed)
            @test occursin("not enabled by default", lower_squashed) ||
                  occursin("does not optimize by default", lower_squashed) ||
                  occursin("does not enable optimization by default", lower_squashed)
        end

        if mentions_issue188 && (
            occursin("global minimum", lower_squashed) ||
            occursin("globally minimal", lower_squashed) ||
            occursin("global optimal", lower_squashed)
        )
            @test occursin("does not claim", lower_squashed) ||
                  occursin("not claim", lower_squashed) ||
                  occursin("no claim", lower_squashed)
        end

        if mentions_issue188 &&
           (occursin("laurent", lower_squashed) || occursin("toricbuilder", lower_squashed)) &&
           (occursin("support", lower_squashed) || occursin("mainline", lower_squashed))
            @test occursin("does not add", lower_squashed) ||
                  occursin("not add", lower_squashed) ||
                  occursin("does not claim", lower_squashed) ||
                  occursin("not claim", lower_squashed) ||
                  occursin("separate", lower_squashed) ||
                  occursin("out of scope", lower_squashed)
        end
    end
end

function _assert_issue188_optimizer_contract(text)
    squashed = _squash_whitespace(text)
    @test occursin(
        "The optional Steinberg factor-count optimizer (#188) is available only through `optimize_elementary_factor_sequence(factors; rules = :safe)`",
        squashed,
    )
    @test occursin("It is not enabled by default", squashed)
    @test occursin(
        "every optimized sequence is accepted only through exact product verification by `verify_steinberg_optimization_certificate`",
        squashed,
    )
    for rule_name in (
        ":identity_removal",
        ":same_position_merge",
        ":inverse_cancellation",
        ":commutator_forward",
        ":commutator_reverse",
        ":disjoint_commutator_identity",
    )
        @test occursin(rule_name, squashed)
    end
    @test occursin(
        "#188 does not change the correctness contract of `elementary_factorization(A)`",
        squashed,
    )
    @test occursin("does not claim global minimum factor counts", squashed)
    @test occursin("does not add Laurent `GL_n` or ToricBuilder support", squashed)
    _assert_issue188_no_overclaims(text)
end

function _assert_issue188_acceptance_audit(text)
    squashed = _squash_whitespace(text)
    for issue_id in ("#288", "#289", "#290", "#291", "#292", "#293")
        @test occursin(issue_id, squashed)
    end
    @test occursin("Park-Woodburn Section 6", squashed)
    @test occursin("| Metric | Original | Optimized | Delta |", text)
    @test occursin("| Factor count | 4 | 1 | -3 |", text)
    @test occursin("| Max monomial degree | 1 | 2 | 1 |", text)
    @test occursin("| Total off-diagonal monomial count | 6 | 2 | -4 |", text)
    @test occursin("| Applied rewrites | `:commutator_forward` | `:commutator_forward` | accepted safe rewrite |", text)
    @test occursin("products_equal = true", squashed)
    @test occursin("verification_status = true", squashed)
    _assert_issue188_optimizer_contract(text)
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

        audit = _read_issue187_acceptance_audit()
        _assert_issue187_acceptance_audit(audit)
        issue188_audit = _read_issue188_acceptance_audit()
        _assert_issue188_acceptance_audit(issue188_audit)
    end

    @testset "ordinary-polynomial Park-Woodburn public contract" begin
        _assert_issue187_public_contract(_read_repo_text(README_PATH))
        _assert_issue187_public_contract(_read_repo_text(DOCS_INDEX_PATH))
        _assert_issue188_optimizer_contract(_read_repo_text(README_PATH))
        _assert_issue188_optimizer_contract(_read_repo_text(DOCS_INDEX_PATH))
    end
end
