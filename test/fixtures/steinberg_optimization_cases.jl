module SteinbergOptimizationFixtureCatalog

using Oscar
using Suslin

const STEINBERG_SECTION_6_REF = "refs/arXiv-alg-geom9405003v1 Section 6"

function _ring_metadata(description, R, generator_names, generators)
    return (;
        description = description,
        object = R,
        generator_names = generator_names,
        generators = generators,
    )
end

function _ordinary_ring_constructor(coefficient, variables)
    return (;
        function_name = :polynomial_ring,
        coefficient = coefficient,
        variables = variables,
    )
end

function _factor_metadata(row, col, coefficient)
    return (;
        row = row,
        col = col,
        coefficient = coefficient,
    )
end

function _product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case(;
    id,
    rule_name,
    description,
    ring_constructor,
    ring,
    matrix_size,
    factor_metadata,
    factors,
    expected_rewrite_factors,
    original_product,
    rewritten_product,
    rewrite_span,
    rule_metadata,
    source_refs,
    consumer_issue_ids,
)
    return (;
        id = id,
        rule_name = rule_name,
        description = description,
        ring_constructor = ring_constructor,
        ring = ring,
        matrix_size = matrix_size,
        factor_metadata = factor_metadata,
        factors = factors,
        expected_rewrite_factors = expected_rewrite_factors,
        original_product = original_product,
        rewritten_product = rewritten_product,
        rewrite_span = rewrite_span,
        rule_metadata = rule_metadata,
        source_refs = source_refs,
        consumer_issue_ids = consumer_issue_ids,
    )
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (;
        id = id,
        base_case_id = base_case_id,
        reason = reason,
    ))
end

function catalog()
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    R_alt, (u, v) = Oscar.polynomial_ring(QQ, ["u", "v"])

    ring_constructor = _ordinary_ring_constructor("QQ", ("x", "y"))
    ring_metadata = _ring_metadata("QQ[x, y]", R, ("x", "y"), (x, y))

    identity_removal_factors = (
        elementary_matrix(3, 1, 2, zero(R), R),
        elementary_matrix(3, 2, 3, x + one(R), R),
    )
    identity_removal_rewrite = (
        elementary_matrix(3, 2, 3, x + one(R), R),
    )
    identity_removal_case = _case(
        id = "steinberg-identity-removal-qq",
        rule_name = :identity_removal,
        description = "Remove an identity elementary factor from an ordinary-polynomial Steinberg rewrite window.",
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        matrix_size = 3,
        factor_metadata = (;
            factors = (
                _factor_metadata(1, 2, zero(R)),
                _factor_metadata(2, 3, x + one(R)),
            ),
            expected_rewrite_factors = (
                _factor_metadata(2, 3, x + one(R)),
            ),
        ),
        factors = identity_removal_factors,
        expected_rewrite_factors = identity_removal_rewrite,
        original_product = _product(identity_removal_factors, R, 3),
        rewritten_product = _product(identity_removal_rewrite, R, 3),
        rewrite_span = (; start = 1, stop = 2),
        rule_metadata = (; indices = (; i = 1, j = 2)),
        source_refs = (
            STEINBERG_SECTION_6_REF,
            "Issue #188 identity-removal fixture",
        ),
        consumer_issue_ids = ("#188",),
    )

    same_position_merge_factors = (
        elementary_matrix(3, 1, 2, x, R),
        elementary_matrix(3, 1, 2, y + one(R), R),
    )
    same_position_merge_rewrite = (
        elementary_matrix(3, 1, 2, x + y + one(R), R),
    )
    same_position_merge_case = _case(
        id = "steinberg-same-position-merge-qq",
        rule_name = :same_position_merge,
        description = "Merge two adjacent elementary factors with the same row and column indices.",
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        matrix_size = 3,
        factor_metadata = (;
            factors = (
                _factor_metadata(1, 2, x),
                _factor_metadata(1, 2, y + one(R)),
            ),
            expected_rewrite_factors = (
                _factor_metadata(1, 2, x + y + one(R)),
            ),
        ),
        factors = same_position_merge_factors,
        expected_rewrite_factors = same_position_merge_rewrite,
        original_product = _product(same_position_merge_factors, R, 3),
        rewritten_product = _product(same_position_merge_rewrite, R, 3),
        rewrite_span = (; start = 1, stop = 2),
        rule_metadata = (; indices = (; i = 1, j = 2)),
        source_refs = (
            STEINBERG_SECTION_6_REF,
            "Issue #188 same-position merge fixture",
        ),
        consumer_issue_ids = ("#188",),
    )

    inverse_cancellation_factors = (
        elementary_matrix(3, 2, 3, x * y + one(R), R),
        elementary_matrix(3, 2, 3, -(x * y + one(R)), R),
    )
    inverse_cancellation_case = _case(
        id = "steinberg-inverse-cancellation-qq",
        rule_name = :inverse_cancellation,
        description = "Cancel consecutive inverse elementary factors with matching indices.",
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        matrix_size = 3,
        factor_metadata = (;
            factors = (
                _factor_metadata(2, 3, x * y + one(R)),
                _factor_metadata(2, 3, -(x * y + one(R))),
            ),
            expected_rewrite_factors = (),
        ),
        factors = inverse_cancellation_factors,
        expected_rewrite_factors = (),
        original_product = _product(inverse_cancellation_factors, R, 3),
        rewritten_product = _product((), R, 3),
        rewrite_span = (; start = 1, stop = 2),
        rule_metadata = (; indices = (; i = 2, j = 3)),
        source_refs = (
            STEINBERG_SECTION_6_REF,
            "Issue #188 inverse-cancellation fixture",
        ),
        consumer_issue_ids = ("#188",),
    )

    commutator_forward_factors = (
        elementary_matrix(3, 1, 2, x + one(R), R),
        elementary_matrix(3, 2, 3, y, R),
        elementary_matrix(3, 1, 2, -(x + one(R)), R),
        elementary_matrix(3, 2, 3, -y, R),
    )
    commutator_forward_rewrite = (
        elementary_matrix(3, 1, 3, (x + one(R)) * y, R),
    )
    commutator_forward_case = _case(
        id = "steinberg-commutator-forward-qq",
        rule_name = :commutator_forward,
        description = "Replay the forward Steinberg commutator relation E_ij(a)E_jl(b)E_ij(-a)E_jl(-b)=E_il(ab).",
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        matrix_size = 3,
        factor_metadata = (;
            factors = (
                _factor_metadata(1, 2, x + one(R)),
                _factor_metadata(2, 3, y),
                _factor_metadata(1, 2, -(x + one(R))),
                _factor_metadata(2, 3, -y),
            ),
            expected_rewrite_factors = (
                _factor_metadata(1, 3, (x + one(R)) * y),
            ),
        ),
        factors = commutator_forward_factors,
        expected_rewrite_factors = commutator_forward_rewrite,
        original_product = _product(commutator_forward_factors, R, 3),
        rewritten_product = _product(commutator_forward_rewrite, R, 3),
        rewrite_span = (; start = 1, stop = 4),
        rule_metadata = (; indices = (; i = 1, j = 2, l = 3)),
        source_refs = (
            STEINBERG_SECTION_6_REF,
            "Issue #188 forward commutator fixture",
        ),
        consumer_issue_ids = ("#188",),
    )

    commutator_reverse_factors = (
        elementary_matrix(3, 2, 3, x, R),
        elementary_matrix(3, 1, 2, y + one(R), R),
        elementary_matrix(3, 2, 3, -x, R),
        elementary_matrix(3, 1, 2, -(y + one(R)), R),
    )
    commutator_reverse_rewrite = (
        elementary_matrix(3, 1, 3, -x * (y + one(R)), R),
    )
    commutator_reverse_case = _case(
        id = "steinberg-commutator-reverse-qq",
        rule_name = :commutator_reverse,
        description = "Replay the reverse Steinberg commutator relation E_ij(a)E_li(b)E_ij(-a)E_li(-b)=E_lj(-ab).",
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        matrix_size = 3,
        factor_metadata = (;
            factors = (
                _factor_metadata(2, 3, x),
                _factor_metadata(1, 2, y + one(R)),
                _factor_metadata(2, 3, -x),
                _factor_metadata(1, 2, -(y + one(R))),
            ),
            expected_rewrite_factors = (
                _factor_metadata(1, 3, -x * (y + one(R))),
            ),
        ),
        factors = commutator_reverse_factors,
        expected_rewrite_factors = commutator_reverse_rewrite,
        original_product = _product(commutator_reverse_factors, R, 3),
        rewritten_product = _product(commutator_reverse_rewrite, R, 3),
        rewrite_span = (; start = 1, stop = 4),
        rule_metadata = (; indices = (; l = 1, i = 2, j = 3)),
        source_refs = (
            STEINBERG_SECTION_6_REF,
            "Issue #188 reverse commutator fixture",
        ),
        consumer_issue_ids = ("#188",),
    )

    disjoint_commutator_identity_factors = (
        elementary_matrix(4, 1, 2, x, R),
        elementary_matrix(4, 3, 4, y + one(R), R),
        elementary_matrix(4, 1, 2, -x, R),
        elementary_matrix(4, 3, 4, -(y + one(R)), R),
    )
    disjoint_commutator_identity_case = _case(
        id = "steinberg-disjoint-commutator-identity-qq",
        rule_name = :disjoint_commutator_identity,
        description = "Record a disjoint commutator window whose Steinberg rewrite is the identity.",
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        matrix_size = 4,
        factor_metadata = (;
            factors = (
                _factor_metadata(1, 2, x),
                _factor_metadata(3, 4, y + one(R)),
                _factor_metadata(1, 2, -x),
                _factor_metadata(3, 4, -(y + one(R))),
            ),
            expected_rewrite_factors = (),
        ),
        factors = disjoint_commutator_identity_factors,
        expected_rewrite_factors = (),
        original_product = _product(disjoint_commutator_identity_factors, R, 4),
        rewritten_product = _product((), R, 4),
        rewrite_span = (; start = 1, stop = 4),
        rule_metadata = (; indices = (; i = 1, j = 2, l = 3, p = 4)),
        source_refs = (
            STEINBERG_SECTION_6_REF,
            "Issue #188 disjoint commutator identity fixture",
        ),
        consumer_issue_ids = ("#188",),
    )

    mismatched_ring_control = _negative_control(
        "steinberg-negative-mismatched-factor-rings",
        same_position_merge_case.id,
        "replace one factor with a matrix over a separate ordinary polynomial ring",
        merge(same_position_merge_case, (;
            factor_metadata = merge(same_position_merge_case.factor_metadata, (;
                factors = (
                    same_position_merge_case.factor_metadata.factors[1],
                    _factor_metadata(1, 2, u + one(R_alt)),
                ),
            )),
            factors = (
                same_position_merge_case.factors[1],
                elementary_matrix(3, 1, 2, u + one(R_alt), R_alt),
            ),
        )),
    )

    stale_product_control = _negative_control(
        "steinberg-negative-stale-expected-product",
        commutator_forward_case.id,
        "store a stale original_product that no longer matches the factor product",
        merge(commutator_forward_case, (;
            original_product = identity_matrix(R, 3),
        )),
    )

    invalid_commutator_indices_control = _negative_control(
        "steinberg-negative-invalid-commutator-indices",
        commutator_forward_case.id,
        "record commutator indices that violate the Section 6 Steinberg relation",
        merge(commutator_forward_case, (;
            rule_metadata = (; indices = (; i = 1, j = 2, l = 1)),
        )),
    )

    return (;
        cases = [
            identity_removal_case,
            same_position_merge_case,
            inverse_cancellation_case,
            commutator_forward_case,
            commutator_reverse_case,
            disjoint_commutator_identity_case,
        ],
        negative_controls = [
            mismatched_ring_control,
            stale_product_control,
            invalid_commutator_indices_control,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
