module SL3MurthyGuptaFixtureCatalog

using Oscar
using Suslin

function _ring_metadata(R, X)
    return _ring_metadata("QQ[X]", R, (X,))
end

function _ring_metadata(description, R, generators)
    return (;
        description = description,
        object = R,
        generator_names = tuple((string(generator) for generator in generators)...),
        generators = generators,
    )
end

function _ring_constructor_metadata(coefficient = "QQ", variables = ("X",))
    return (;
        function_name = :polynomial_ring,
        coefficient = coefficient,
        variables = variables,
    )
end

function _target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _case(id, branch, variable, entries, target, witness, expected_current_solver, consumer_issue_ids;
        ring_constructor = _ring_constructor_metadata(),
        ring = _ring_metadata(parent(variable), variable),
        extra = (;))
    base = (;
        id = id,
        branch = branch,
        ring_constructor = ring_constructor,
        ring = ring,
        variable = variable,
        entries = entries,
        target = target,
        murthy_path = true,
        expected_current_solver = expected_current_solver,
        witnesses = witness,
        source_refs = ("Park-Woodburn arXiv:alg-geom/9405003 section 5",),
        consumer_issue_ids = consumer_issue_ids,
    )
    return merge(base, extra)
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (; id, base_case_id, reason))
end

function catalog()
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    now_supported = (; status = :passes)
    expected_local_failure = (; status = :staged_fail, message_substring = "staged local SL_3 solver failure")
    local_consumer_issues = ("#182", "#208", "#207", "#209", "#210")

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
        now_supported,
        ("#74",),
    )

    q0_nonunit_extracted_bezout_case = _case(
        "mg-q0-nonunit-extracted-bezout-resultant",
        :q0_nonunit_bezout_resultant,
        X,
        (;
            p = X^3 + X + 1,
            q = X^2,
            r = X^3 - X^2 + 2 * X,
            s = X^2 - X + 1,
        ),
        _target(R, X^3 + X + 1, X^2, X^3 - X^2 + 2 * X, X^2 - X + 1),
        ((
            p0 = one(R),
            q0 = zero(R),
            p_prime = one(R) - X,
            q_prime = -X^2 + X - 1,
            resultant = one(R),
            p_prime_degree = 1,
            q_prime_degree = 2,
            branch_unit = one(R),
            case1_entries = (;
                p = X^3 - X^2 + 2 * X,
                q = X^2 - X + 1,
                r = -X^2 + X - 1,
                s = one(R) - X,
            ),
        ),),
        now_supported,
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

    RU, (u, UX) = Oscar.polynomial_ring(QQ, ["u", "X"])
    local_ring = _ring_metadata("QQ[u, X]", RU, (u, UX))
    local_ring_constructor = _ring_constructor_metadata("QQ", ("u", "X"))
    local_context = (;
        kind = :localization_at_maximal_ideal,
        description = "QQ[u] localized at (u), with X as the Section 5 variable",
        selected_variable = UX,
        maximal_ideal_generators = (u,),
        residue_description = "u => 0",
    )

    function _local_unit(unit, residue_unit, residue_inverse, generators, coefficients)
        return (;
            context = merge(local_context, (; maximal_ideal_generators = generators)),
            unit,
            residue_unit,
            residue_inverse,
            maximal_ideal_generators = generators,
            residue_difference_coefficients = coefficients,
            global_unit = false,
        )
    end

    local_q_degree_entries = (; p = UX^2 + u * UX + one(RU), q = UX * (UX^2 + u * UX + one(RU)) + one(RU), r = -one(RU), s = -UX)
    local_q_degree_case = _case(
        "mg-local-q-degree-qq-u-x",
        :q_degree_normalization,
        UX,
        local_q_degree_entries,
        _target(RU, local_q_degree_entries.p, local_q_degree_entries.q, local_q_degree_entries.r, local_q_degree_entries.s),
        ((
            quotient = UX,
            remainder = one(RU),
            normalized_s = zero(RU),
        ),),
        expected_local_failure,
        local_consumer_issues;
        ring_constructor = local_ring_constructor,
        ring = local_ring,
        extra = (;
            local_contract = true,
            requires_local_units = false,
            requires_bezout_witness = false,
        ),
    )

    local_q0_unit_entries = (; p = UX^2 + (u + one(RU)) * UX + one(RU), q = UX + u + one(RU), r = UX + (UX^2 + (u + one(RU)) * UX + one(RU)) * UX, s = UX^2 + (u + one(RU)) * UX + one(RU))
    local_q0_unit_case = _case(
        "mg-local-q0-unit-at-u",
        :q0_unit_recursion,
        UX,
        local_q0_unit_entries,
        _target(RU, local_q0_unit_entries.p, local_q0_unit_entries.q, local_q0_unit_entries.r, local_q0_unit_entries.s),
        ((
            p0 = one(RU),
            q0 = u + one(RU),
            local_unit_witness = _local_unit(u + one(RU), one(RU), one(RU), (u,), (one(RU),)),
            formal_right_e21_coefficient = "-1/(1 + u)",
        ),),
        expected_local_failure,
        local_consumer_issues;
        ring_constructor = local_ring_constructor,
        ring = local_ring,
        extra = (;
            local_contract = true,
            requires_local_units = true,
            requires_bezout_witness = false,
        ),
    )

    local_q0_nonunit_entries = (; p = UX^2 + u * UX + one(RU), q = UX + u, r = UX + (UX^2 + u * UX + one(RU)) * UX, s = UX^2 + u * UX + one(RU))
    local_q0_nonunit_case = _case(
        "mg-local-q0-nonunit-bezout-at-u",
        :q0_nonunit_bezout_resultant,
        UX,
        local_q0_nonunit_entries,
        _target(RU, local_q0_nonunit_entries.p, local_q0_nonunit_entries.q, local_q0_nonunit_entries.r, local_q0_nonunit_entries.s),
        ((
            p0 = one(RU),
            q0 = u,
            p_prime = one(RU),
            q_prime = UX,
            resultant = one(RU),
            p_prime_degree = 0,
            q_prime_degree = 1,
            branch_unit = u + one(RU),
            branch_unit_witness = _local_unit(u + one(RU), one(RU), one(RU), (u,), (one(RU),)),
            case1_entries = (;
                p = local_q0_nonunit_entries.p + UX,
                q = local_q0_nonunit_entries.q + one(RU),
                r = UX,
                s = one(RU),
            ),
        ),),
        expected_local_failure,
        local_consumer_issues;
        ring_constructor = local_ring_constructor,
        ring = local_ring,
        extra = (;
            local_contract = true,
            requires_local_units = true,
            requires_bezout_witness = true,
        ),
    )

    nonmonic_entries = (; p = 2 * X + one(R), q = X, r = R(2), s = one(R))
    nonmonic = _case(
        "mg-negative-nonmonic-p",
        :open_slice_control,
        X,
        nonmonic_entries,
        _target(R, nonmonic_entries.p, nonmonic_entries.q, nonmonic_entries.r, nonmonic_entries.s),
        (),
        now_supported,
        ("#206",),
    )

    det_bad_entries = (; p = X + one(R), q = zero(R), r = zero(R), s = one(R))
    det_bad = _case(
        "mg-negative-determinant-not-one",
        :open_slice_control,
        X,
        det_bad_entries,
        _target(R, det_bad_entries.p, det_bad_entries.q, det_bad_entries.r, det_bad_entries.s),
        (),
        now_supported,
        ("#206",),
    )

    split_negative = _negative_control(
        "mg-negative-corrupted-split-witness",
        split_lemma_case.id,
        "split witness a no longer reconstructs the target",
        merge(split_lemma_case, (; witnesses = (merge(first(split_lemma_case.witnesses), (; a = first(split_lemma_case.witnesses).a + one(R))),))),
    )

    local_unit_negative = _negative_control(
        "mg-negative-corrupted-local-unit-witness",
        local_q0_unit_case.id,
        "local unit residue equation is corrupted",
        merge(local_q0_unit_case, (; witnesses = (merge(first(local_q0_unit_case.witnesses), (;
            local_unit_witness = merge(first(local_q0_unit_case.witnesses).local_unit_witness, (;
                residue_difference_coefficients = (zero(RU),),
            )),
        )),))),
    )

    bezout_negative = _negative_control(
        "mg-negative-corrupted-bezout-equality",
        local_q0_nonunit_case.id,
        "p_prime*p - q_prime*q no longer equals the resultant",
        merge(local_q0_nonunit_case, (; witnesses = (merge(first(local_q0_nonunit_case.witnesses), (; p_prime = one(RU) + u)),))),
    )

    return (;
        ring = _ring_metadata(R, X),
        cases = [
            q_degree_normalization_case,
            split_lemma_case,
            q0_unit_case,
            q0_nonunit_normalizes_case,
            q0_nonunit_normalized_bezout_case,
            q0_nonunit_extracted_bezout_case,
            open_slice_case,
            local_q_degree_case,
            local_q0_unit_case,
            local_q0_nonunit_case,
        ],
        negative_controls = [
            nonmonic,
            det_bad,
            split_negative,
            local_unit_negative,
            bezout_negative,
        ],
    )
end

end
