module ToricBuilderIssue38Cases

using Oscar
using Suslin

const ISSUE38_FAILURE_SUBSTRINGS = (
    "staged SL_n to local SL_3 reduction failure",
    "failed to solve local SL_3 obligation",
)

function _issue38_q_block(R, u, v)
    return matrix(R, [
        1 + v^-1              1      1  0  1 + v^-1              1 + v^-1;
        u*v^-1 + 1 + v^-1     u      1  1  u*v^-1 + 1 + v^-1     u*v^-1 + v^-1;
        u*v^-1                u*v    0  0  u*v^-1                u*v^-1;
        1                     1      1  0  1                     1;
        1 + v^-1              v      0  0  v^-1                  1 + v^-1;
        u + v                 u*v    v  v  u + v                 u + v
    ])
end

_ring_metadata(R, u, v) = (;
    description = "GF(2)[u^+/-1, v^+/-1]",
    object = R,
    generators = (u, v),
    variables = ("u", "v"),
)

_issue38_provenance() = (;
    source = :toricbuilder_issue_38_mwe,
    issue = "#38",
    issue_url = "https://github.com/nzy1997/Suslin.jl/issues/38",
    source_description = "upper-left Q block of transfer_result.column_transformation for the 1 + x + x*y color-code example",
    reported_main_commit = "c985b1aac9fc9152d860e4e90d012964867bb27d",
)

_expected_failure_status() = (;
    kind = :staged_sl3_local_obligation_failure,
    message_substrings = ISSUE38_FAILURE_SUBSTRINGS,
)

function catalog()
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    Q = _issue38_q_block(R, u, v)

    determinant_profile = Suslin.classify_laurent_determinant(Q)
    row_normalization = Suslin.normalize_laurent_gl_matrix(Q)
    Dcol = diagonal_matrix(R, [inv(determinant_profile.determinant), one(R), one(R), one(R), one(R), one(R)])
    col_core = Q * Dcol

    row_status = _expected_failure_status()
    column_status = _expected_failure_status()

    return (;
        cases = [
            (;
                id = "toricbuilder-issue-38-q-block",
                kind = :toricbuilder_issue38_q_block,
                ring = _ring_metadata(R, u, v),
                dimensions = (; matrix = (6, 6)),
                inputs = (; matrix = Q),
                determinant_profile = (;
                    expected_determinant = determinant_profile.determinant,
                    expected_class = determinant_profile.classification,
                    monomial_exponents = determinant_profile.monomial_exponents,
                    monomial_coefficient = determinant_profile.monomial_coefficient,
                ),
                normalizations = (;
                    row = (;
                        core = row_normalization.normalized_matrix,
                        normalization = row_normalization,
                        expected_current_status = row_status,
                    ),
                    column = (;
                        core = col_core,
                        correction_factor = Dcol,
                        expected_current_status = column_status,
                    ),
                ),
                expected_current_status = (;
                    row = row_status,
                    column = column_status,
                ),
                provenance = _issue38_provenance(),
                consumer_test_ids = ("issue-39-toricbuilder-issue38-fixture",),
            ),
        ],
    )
end

end
