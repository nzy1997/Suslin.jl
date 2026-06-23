module QuillenPatchFixtureCatalog

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
    kind,
    stage_coverage,
    ring_constructor,
    ring,
    size,
    substitution_variable,
    target_matrix,
    base_matrix,
    denominator_data,
    local_factors,
    expected,
    patched_substitution_witness,
    source_refs,
    consumer_issue_ids,
)
    return (;
        id = id,
        kind = kind,
        stage_coverage = stage_coverage,
        ring_constructor = ring_constructor,
        ring = ring,
        size = size,
        substitution_variable = substitution_variable,
        target_matrix = target_matrix,
        base_matrix = base_matrix,
        denominator_data = denominator_data,
        local_factors = local_factors,
        expected = expected,
        patched_substitution_witness = patched_substitution_witness,
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

function _denominator_data(denominator, coverage_multiplier)
    return (;
        denominator = denominator,
        coverage_multiplier = coverage_multiplier,
    )
end

function _certificate(indices, denominators)
    return (;
        indices = indices,
        denominators = denominators,
    )
end

function _correction(row, col, entry)
    return (;
        row = row,
        col = col,
        entry = entry,
    )
end

function _local_factor(;
    certificate,
    denominator,
    coverage_multiplier,
    correction,
    factor,
    expected_correction,
)
    return (;
        certificate = certificate,
        denominator = denominator,
        coverage_multiplier = coverage_multiplier,
        correction = correction,
        factor = factor,
        expected_correction = expected_correction,
    )
end

function _elementary_local_factor(n, row, col, entry, denominator, coverage_multiplier, R)
    weighted_entry = coverage_multiplier * denominator * entry
    expected_correction = elementary_matrix(n, row, col, weighted_entry, R)
    return _local_factor(
        certificate = _certificate([row, col], [denominator, denominator]),
        denominator = denominator,
        coverage_multiplier = coverage_multiplier,
        correction = _correction(row, col, entry),
        factor = expected_correction,
        expected_correction = expected_correction,
    )
end

function _product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _expected_pass(global_correction)
    return (;
        current_status = :passes,
        global_correction = global_correction,
    )
end

function catalog()
    RQ, (X, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    R2, (X2, r2, s2) = Oscar.polynomial_ring(GF(2), ["X", "r", "s"])

    qq_ring_constructor = _ordinary_ring_constructor("QQ", ("X", "r", "g"))
    gf2_ring_constructor = _ordinary_ring_constructor("GF(2)", ("X", "r", "s"))

    qq_ring = _ring_metadata("QQ[X, r, g]", RQ, ("X", "r", "g"), (X, r, g))
    gf2_ring = _ring_metadata("GF(2)[X, r, s]", R2, ("X", "r", "s"), (X2, r2, s2))

    n = 3

    two_open_entry = X + g + one(RQ)
    two_open_factors = (
        _elementary_local_factor(n, 1, 2, two_open_entry, r, one(RQ), RQ),
        _elementary_local_factor(n, 1, 2, two_open_entry, one(RQ) - r, one(RQ), RQ),
    )
    two_open_target = _product([factor.factor for factor in two_open_factors], RQ, n)
    two_open_case = _case(
        id = "quillen-two-open-cover-qq",
        kind = :two_open_cover,
        stage_coverage = :supported,
        ring_constructor = qq_ring_constructor,
        ring = qq_ring,
        size = n,
        substitution_variable = X,
        target_matrix = two_open_target,
        base_matrix = matrix(RQ, [
            one(RQ)  X        g;
            zero(RQ) one(RQ)  r;
            zero(RQ) zero(RQ) one(RQ)
        ]),
        denominator_data = (
            _denominator_data(r, one(RQ)),
            _denominator_data(one(RQ) - r, one(RQ)),
        ),
        local_factors = two_open_factors,
        expected = _expected_pass(two_open_target),
        patched_substitution_witness = nothing,
        source_refs = ("test/expert/quillen_patching_exact.jl ordinary QQ two-open example",),
        consumer_issue_ids = ("#63", "#99"),
    )

    nontrivial_entry = X + r + g
    first_multiplier = one(RQ) + g * (one(RQ) - r)
    second_multiplier = one(RQ) - g * r
    nontrivial_factors = (
        _elementary_local_factor(n, 1, 3, nontrivial_entry, r, first_multiplier, RQ),
        _elementary_local_factor(n, 1, 3, nontrivial_entry, one(RQ) - r, second_multiplier, RQ),
    )
    nontrivial_target = _product([factor.factor for factor in nontrivial_factors], RQ, n)
    nontrivial_case = _case(
        id = "quillen-nontrivial-multipliers-qq",
        kind = :nontrivial_coverage_multipliers,
        stage_coverage = :supported,
        ring_constructor = qq_ring_constructor,
        ring = qq_ring,
        size = n,
        substitution_variable = X,
        target_matrix = nontrivial_target,
        base_matrix = identity_matrix(RQ, n),
        denominator_data = (
            _denominator_data(r, first_multiplier),
            _denominator_data(one(RQ) - r, second_multiplier),
        ),
        local_factors = nontrivial_factors,
        expected = _expected_pass(nontrivial_target),
        patched_substitution_witness = nothing,
        source_refs = ("Issue 99 nontrivial coverage multiplier acceptance case",),
        consumer_issue_ids = ("#63", "#99"),
    )

    supplied_entry = X2 + s2 + one(R2)
    supplied_factors = (
        _elementary_local_factor(n, 1, 3, supplied_entry, r2, one(R2), R2),
        _elementary_local_factor(n, 1, 3, supplied_entry, one(R2) + r2, one(R2), R2),
    )
    supplied_target = _product([factor.factor for factor in supplied_factors], R2, n)
    supplied_certificate_case = _case(
        id = "quillen-supplied-local-certificate-gf2",
        kind = :supplied_local_certificate,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = gf2_ring,
        size = n,
        substitution_variable = X2,
        target_matrix = supplied_target,
        base_matrix = identity_matrix(R2, n),
        denominator_data = (
            _denominator_data(r2, one(R2)),
            _denominator_data(one(R2) + r2, one(R2)),
        ),
        local_factors = supplied_factors,
        expected = _expected_pass(supplied_target),
        patched_substitution_witness = nothing,
        source_refs = ("Issue 99 GF(2) supplied-local-certificate acceptance case",),
        consumer_issue_ids = ("#63", "#99"),
    )

    shifted_X = X + r^2 * g
    patched_source = elementary_matrix(n, 1, 2, X + g + one(RQ), RQ)
    patched_expected = elementary_matrix(n, 1, 2, shifted_X + g + one(RQ), RQ)
    patched_entry = shifted_X + g + one(RQ)
    patched_factors = (
        _elementary_local_factor(n, 1, 2, patched_entry, r, one(RQ), RQ),
        _elementary_local_factor(n, 1, 2, patched_entry, one(RQ) - r, one(RQ), RQ),
    )
    patched_target = _product([factor.factor for factor in patched_factors], RQ, n)
    patched_witness_case = _case(
        id = "quillen-patched-substitution-witness-qq",
        kind = :patched_substitution_witness,
        stage_coverage = :supported,
        ring_constructor = qq_ring_constructor,
        ring = qq_ring,
        size = n,
        substitution_variable = X,
        target_matrix = patched_target,
        base_matrix = patched_source,
        denominator_data = (
            _denominator_data(r, one(RQ)),
            _denominator_data(one(RQ) - r, one(RQ)),
        ),
        local_factors = patched_factors,
        expected = _expected_pass(patched_target),
        patched_substitution_witness = (;
            matrix = patched_source,
            variable = X,
            denominator = r,
            exponent = 2,
            shift = g,
            expected_matrix = patched_expected,
        ),
        source_refs = ("test/expert/quillen_induction.jl patched_substitution witness shape",),
        consumer_issue_ids = ("#63", "#99", "#102"),
    )

    constructive_entry = X2 * r2 + s2 + one(R2)
    constructive_factors = (
        _elementary_local_factor(n, 2, 1, constructive_entry, s2, one(R2), R2),
        _elementary_local_factor(n, 2, 1, constructive_entry, one(R2) + s2, one(R2), R2),
    )
    constructive_target = _product([factor.factor for factor in constructive_factors], R2, n)
    constructive_case = _case(
        id = "quillen-constructive-acceptance-gf2",
        kind = :constructive_acceptance,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = gf2_ring,
        size = n,
        substitution_variable = X2,
        target_matrix = constructive_target,
        base_matrix = identity_matrix(R2, n),
        denominator_data = (
            _denominator_data(s2, one(R2)),
            _denominator_data(one(R2) + s2, one(R2)),
        ),
        local_factors = constructive_factors,
        expected = _expected_pass(constructive_target),
        patched_substitution_witness = nothing,
        source_refs = ("Issue 99 GF(2) constructive Quillen acceptance case",),
        consumer_issue_ids = ("#63", "#99"),
    )

    bad_cover = _negative_control(
        "quillen-uncovered-denominator-control",
        "quillen-two-open-cover-qq",
        "mutated coverage multiplier",
        merge(two_open_case, (;
            denominator_data = (
                two_open_case.denominator_data[1],
                merge(two_open_case.denominator_data[2], (; coverage_multiplier = r)),
            ),
        )),
    )

    bad_factor = _negative_control(
        "quillen-tampered-local-factor-control",
        "quillen-supplied-local-certificate-gf2",
        "mutated local factor matrix",
        merge(supplied_certificate_case, (;
            local_factors = (
                merge(supplied_certificate_case.local_factors[1], (;
                    factor = supplied_certificate_case.local_factors[1].factor *
                        elementary_matrix(supplied_certificate_case.size, 1, 3, one(R2), R2),
                )),
                supplied_certificate_case.local_factors[2],
            ),
        )),
    )

    return (;
        cases = [
            two_open_case,
            nontrivial_case,
            supplied_certificate_case,
            patched_witness_case,
            constructive_case,
        ],
        negative_controls = [
            bad_cover,
            bad_factor,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
