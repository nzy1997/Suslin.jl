module ParkWoodburnPolynomialNormalityFixtureCatalog

using Oscar
using Suslin

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

function _case(;
        id,
        section2_layer,
        ring_constructor,
        ring,
        inputs,
        target_matrix,
        expected_convention,
        source_refs,
        consumer_issue_ids)
    return (;
        id = id,
        section2_layer = section2_layer,
        ring_constructor = ring_constructor,
        ring = ring,
        inputs = inputs,
        target_matrix = target_matrix,
        expected_convention = expected_convention,
        source = "Park-Woodburn issue #190",
        source_refs = source_refs,
        provenance = (;
            source = "Park-Woodburn issue #190",
            issue = "#190",
            section2_layer = section2_layer,
        ),
        consumer_issue_ids = consumer_issue_ids,
    )
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(
        entry,
        (;
            id = id,
            base_case_id = base_case_id,
            reason = reason,
        ),
    )
end

function catalog()
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    ring_constructor = _ordinary_ring_constructor("QQ", ("x", "y"))
    ring_metadata = _ring_metadata("QQ[x, y]", R, ("x", "y"), (x, y))

    cohn_target = identity_matrix(R, 3)
    cohn_i = 1
    cohn_j = 2
    cohn_a = x + y
    cohn_v = [one(R), x, y]
    vi, vj = cohn_v[cohn_i], cohn_v[cohn_j]
    for row in 1:3
        cohn_target[row, cohn_i] += cohn_a * cohn_v[row] * vj
        cohn_target[row, cohn_j] -= cohn_a * cohn_v[row] * vi
    end
    cohn_case = _case(
        id = "pw-section2-cohn-type-qq",
        section2_layer = :cohn_type,
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        inputs = (;
            i = cohn_i,
            j = cohn_j,
            a = cohn_a,
            v = cohn_v,
        ),
        target_matrix = cohn_target,
        expected_convention = :cohn_type,
        source_refs = ("Park-Woodburn issue 190 Cohn-type section 2 fixture",),
        consumer_issue_ids = ("#190",),
    )

    cohn_tampered_target = copy(cohn_target)
    cohn_tampered_target[1, 1] += one(R)
    cohn_control = _negative_control(
        "pw-section2-cohn-type-tampered-target-control",
        cohn_case.id,
        "tampered cohn target does not match i/j/a/v inputs",
        merge(cohn_case, (; target_matrix = cohn_tampered_target)),
    )

    rank_one_v = [one(R), x, y]
    rank_one_w = [-x, one(R), zero(R)]
    rank_one_g = [one(R), zero(R), zero(R)]
    rank_one_target = identity_matrix(R, 3)
    for row in 1:3, col in 1:3
        rank_one_target[row, col] += rank_one_v[row] * rank_one_w[col]
    end
    rank_one_case = _case(
        id = "pw-section2-orthogonal-rank-one-qq",
        section2_layer = :orthogonal_rank_one,
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        inputs = (;
            v = rank_one_v,
            w = rank_one_w,
            g = rank_one_g,
        ),
        target_matrix = rank_one_target,
        expected_convention = :orthogonal_rank_one,
        source_refs = ("Park-Woodburn issue 190 rank-one section 2 fixture",),
        consumer_issue_ids = ("#190",),
    )

    rank_one_bad_control_case = _negative_control(
        "pw-section2-rank-one-bad-orthogonality-control",
        rank_one_case.id,
        "rank-one control with non-orthogonal v and w",
        merge(rank_one_case, (; inputs = (;
            v = rank_one_v,
            w = [rank_one_w[1], rank_one_w[2], rank_one_w[3] + one(R)],
            g = rank_one_g,
        ))),
    )

    conjugated_B = matrix(R, [one(R) zero(R) zero(R); x one(R) zero(R); zero(R) one(R) one(R)])
    conjugated_i = 1
    conjugated_j = 3
    conjugated_a = x + y + one(R)
    conjugated_E = elementary_matrix(3, conjugated_i, conjugated_j, conjugated_a, R)
    conjugated_target = conjugated_B * conjugated_E * inv(conjugated_B)
    conjugated_case = _case(
        id = "pw-section2-conjugated-elementary-qq",
        section2_layer = :conjugated_elementary,
        ring_constructor = ring_constructor,
        ring = ring_metadata,
        inputs = (;
            B = conjugated_B,
            i = conjugated_i,
            j = conjugated_j,
            a = conjugated_a,
        ),
        target_matrix = conjugated_target,
        expected_convention = :conjugated_elementary,
        source_refs = ("Park-Woodburn issue 190 conjugated elementary section 2 fixture",),
        consumer_issue_ids = ("#190",),
    )

    conjugated_tampered_target = copy(conjugated_target)
    conjugated_tampered_target[2, 2] += one(R)
    conjugated_control = _negative_control(
        "pw-section2-conjugated-elementary-tampered-target-control",
        conjugated_case.id,
        "tampered conjugated target does not match conjugated elementary decomposition",
        merge(conjugated_case, (; target_matrix = conjugated_tampered_target)),
    )

    return (;
        cases = [
            cohn_case,
            rank_one_case,
            conjugated_case,
        ],
        negative_controls = [
            cohn_control,
            rank_one_bad_control_case,
            conjugated_control,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

function negative_controls_by_id()
    return Dict(entry.id => entry for entry in catalog().negative_controls)
end

end
