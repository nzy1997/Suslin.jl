module LaurentLazyDeterminantCases

using Oscar
using Suslin

include("toricbuilder_issue38_cases.jl")

_matrix_dimensions(A) = (nrows(A), ncols(A))

function _ring_metadata(R, generators, variables, description)
    return (;
        description,
        object = R,
        generators,
        variables,
    )
end

_synthetic_provenance(description) = (;
    source = :synthetic,
    issue = "#154",
    description,
)

function _determinant_profile(A)
    profile = Suslin.classify_laurent_determinant(A)
    return (;
        expected_determinant = profile.determinant,
        expected_class = profile.classification,
        monomial_exponents = profile.monomial_exponents,
        monomial_coefficient = profile.monomial_coefficient,
    )
end

function _expected_correction(profile)
    if profile.expected_class == :one
        return (;
            supported = true,
            kind = :identity,
            supports = (:row_core, :column_core),
            unsupported_reason = nothing,
        )
    elseif profile.expected_class == :laurent_monomial_unit
        return (;
            supported = true,
            kind = :monomial_unit_diagonal,
            supports = (:row_core, :column_core),
            unsupported_reason = nothing,
        )
    else
        return (;
            supported = false,
            kind = :unsupported,
            supports = (),
            unsupported_reason = profile.expected_class,
        )
    end
end

function _supported_normalizations(A)
    R = base_ring(A)
    n = nrows(A)
    row = Suslin.normalize_laurent_gl_matrix(A)
    determinant = det(A)
    column_correction = diagonal_matrix(R, [i == 1 ? inv(determinant) : one(R) for i in 1:n])
    return (;
        row = (; core = row.normalized_matrix, normalization = row),
        column = (; core = A * column_correction, correction_factor = column_correction),
    )
end

function _negative_control(base_case_id, expected_failure)
    return (;
        kind = :metadata_mutation,
        base_case_id,
        expected_failure,
    )
end

function _issue38_entry()
    issue38 = only(ToricBuilderIssue38Cases.catalog().cases)
    return (;
        id = "issue-38-q-block-lazy-determinant",
        kind = :issue38_q_block_lazy_determinant,
        ring = issue38.ring,
        dimensions = issue38.dimensions,
        inputs = issue38.inputs,
        determinant_profile = issue38.determinant_profile,
        expected_correction = _expected_correction(issue38.determinant_profile),
        normalizations = issue38.normalizations,
        negative_control = _negative_control(
            "issue-38-q-block-lazy-determinant",
            :determinant_class_metadata_mismatch,
        ),
        provenance = (;
            source = :wrapped_fixture,
            issue = "#154",
            source_fixture_id = issue38.id,
            source_issue = issue38.provenance.issue,
            description = "reuse ToricBuilder Issue #38 Q-block fixture without copying matrix data",
        ),
        consumer_test_ids = ("issue-154-lazy-laurent-determinant-fixtures",),
    )
end

function _synthetic_entry(id, kind, ring, A, description)
    profile = _determinant_profile(A)
    return (;
        id,
        kind,
        ring,
        dimensions = (; matrix = _matrix_dimensions(A)),
        inputs = (; matrix = A),
        determinant_profile = profile,
        expected_correction = _expected_correction(profile),
        normalizations = profile.expected_class in (:one, :laurent_monomial_unit) ?
            _supported_normalizations(A) : nothing,
        negative_control = _negative_control(id, :determinant_class_metadata_mismatch),
        provenance = _synthetic_provenance(description),
        consumer_test_ids = ("issue-154-lazy-laurent-determinant-fixtures",),
    )
end

function catalog()
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    ring = _ring_metadata(R, (x, y), ("x", "y"), "GF(2)[x^+/-1, y^+/-1]")

    determinant_one = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])
    monomial_unit = matrix(R, [
        x^-1 * y one(R);
        zero(R) one(R)
    ])
    non_unit = matrix(R, [
        x + one(R) zero(R);
        zero(R) one(R)
    ])

    return (;
        cases = [
            _issue38_entry(),
            _synthetic_entry(
                "determinant-one-triangular",
                :determinant_one_triangular,
                ring,
                determinant_one,
                "triangular Laurent matrix with determinant one and identity correction",
            ),
            _synthetic_entry(
                "monomial-unit-row-column-cores",
                :monomial_unit_row_column_cores,
                ring,
                monomial_unit,
                "triangular Laurent matrix with monomial-unit determinant and explicit row/column cores",
            ),
            _synthetic_entry(
                "non-unit-determinant-negative",
                :non_unit_determinant_negative,
                ring,
                non_unit,
                "triangular Laurent matrix whose non-unit determinant must be rejected",
            ),
        ],
    )
end

end
