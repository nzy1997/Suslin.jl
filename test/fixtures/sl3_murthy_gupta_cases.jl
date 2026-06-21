module SL3MurthyGuptaFixtureCatalog

using Oscar
using Suslin

function _ring_metadata(R, X)
    return (;
        description = "QQ[X]",
        object = R,
        generators = (X,),
    )
end

function _ring_constructor_metadata()
    return (;
        function_name = :polynomial_ring,
        coefficient = "QQ",
        variables = ("X",),
    )
end

function _target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _case(id, branch, variable, entries, target, witness, expected_current_solver, consumer_issue_ids)
    return (
        id = id,
        branch = branch,
        ring_constructor = _ring_constructor_metadata(),
        ring = _ring_metadata(parent(variable), variable),
        variable = variable,
        entries = entries,
        target = target,
        murthy_path = true,
        expected_current_solver = expected_current_solver,
        witnesses = witness,
        source_refs = ("Park-Woodburn arXiv:alg-geom/9405003 section 5",),
        consumer_issue_ids = consumer_issue_ids,
    )
end

function catalog()
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    common_failure = (; status = :staged_fail, message_substring = "staged local SL_3 solver failure")
    now_supported = (; status = :passes)

    q_degree_normalization_case = _case(
        "mg-q-degree-normalization",
        :q_degree_normalization,
        X,
        (;
            p = X^2 + 1,
            q = X^3 + X + 1,
            r = -one(R),
            s = -X,
        ),
        _target(R, X^2 + 1, X^3 + X + 1, -one(R), -X),
        ((
            quotient = X,
            remainder = one(R),
            normalized_s = zero(R),
        ),),
        now_supported,
        ("#71",),
    )

    split_lemma_case = _case(
        "mg-split-lemma-x-square",
        :split_lemma,
        X,
        (;
            p = X^2,
            q = one(R),
            r = X^3 + X^2 - 1,
            s = X + 1,
        ),
        _target(R, X^2, one(R), X^3 + X^2 - 1, X + 1),
        ((
            a = X,
            a_prime = X,
            b = one(R),
            c = X^3 + X^2 - 1,
            c1 = X^3 + X^2 - 1,
            c2 = X^3 + X^2 - 1,
            d1 = X^2 + X,
            d2 = X^2 + X,
            d = X + 1,
        ),),
        now_supported,
        ("#72",),
    )

    q0_unit_case = _case(
        "mg-q0-unit-recursion",
        :q0_unit_recursion,
        X,
        (;
            p = X + 1,
            q = one(R),
            r = X^2 + 2 * X,
            s = X + 1,
        ),
        _target(R, X + 1, one(R), X^2 + 2 * X, X + 1),
        ((
            p0 = one(R),
            q0 = one(R),
            q0_inverse = one(R),
            right_e21_coefficient = -one(R),
            normalized_p = X,
            normalized_r = X^2 + X - 1,
            normalized_s = X + 1,
            split = (;
                a = X,
                a_prime = one(R),
                b = one(R),
                c = X^2 + X - 1,
                c1 = X^2 + X - 1,
                c2 = X^2 + X - 1,
                d1 = X + 1,
                d2 = X^2 + X,
                d = X + 1,
            ),
        ),),
        now_supported,
        ("#73",),
    )

    q0_nonunit_normalizes_case = _case(
        "mg-q0-nonunit-normalizes-to-q0-unit",
        :q_degree_normalization,
        X,
        (;
            p = X + 1,
            q = X,
            r = X + 2,
            s = X + 1,
        ),
        _target(R, X + 1, X, X + 2, X + 1),
        ((
            quotient = one(R),
            remainder = -one(R),
            normalized_s = -one(R),
        ),),
        now_supported,
        ("#71", "#73"),
    )

    q0_nonunit_normalized_bezout_case = _case(
        "mg-q0-nonunit-normalized-bezout-resultant",
        :q0_nonunit_bezout_resultant,
        X,
        (;
            p = X^2 + 1,
            q = X,
            r = X^2 + X + 1,
            s = X + 1,
        ),
        _target(R, X^2 + 1, X, X^2 + X + 1, X + 1),
        ((
            p0 = one(R),
            q0 = zero(R),
            p_prime = one(R),
            q_prime = X,
            resultant = one(R),
            p_prime_degree = 0,
            q_prime_degree = 1,
            branch_unit = one(R),
            case1_entries = (;
                p = X^2 + X + 1,
                q = X + 1,
                r = X,
                s = one(R),
            ),
        ),),
        common_failure,
        ("#74",),
    )

    open_slice_case = _case(
        "mg-open-slice-control",
        :open_slice_control,
        X,
        (;
            p = X + 1,
            q = one(R),
            r = X,
            s = one(R),
        ),
        _target(R, X + 1, one(R), X, one(R)),
        (),
        (; status = :passes),
        ("#75",),
    )

    return (;
        ring = _ring_metadata(R, X),
        cases = [
            q_degree_normalization_case,
            split_lemma_case,
            q0_unit_case,
            q0_nonunit_normalizes_case,
            q0_nonunit_normalized_bezout_case,
            open_slice_case,
        ],
    )
end

end
