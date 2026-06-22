module ECPColumnFixtureCatalog

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

function _case(; id, kind, stage_coverage, ring_constructor, ring, variable_order, entries, column_order, monicity, witnesses, expected, source_refs, consumer_issue_ids)
    return (;
        id = id,
        kind = kind,
        stage_coverage = stage_coverage,
        ring_constructor = ring_constructor,
        ring = ring,
        variable_order = variable_order,
        entries = entries,
        column_order = column_order,
        monicity = monicity,
        witnesses = witnesses,
        expected = expected,
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

function _identity_monicity(selected_entry, variable_name, transformed_entry)
    return (;
        selected_entry = selected_entry,
        variable_name = variable_name,
        substitution = NamedTuple(),
        transformed_entry = transformed_entry,
    )
end

function _column_values(entry)
    return entry.entries
end

function catalog()
    R2, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    RQ, (X, Y) = Oscar.polynomial_ring(QQ, ["x", "y"])

    gf2_ring_constructor = _ordinary_ring_constructor("GF(2)", ("x", "y"))
    qq_ring_constructor = _ordinary_ring_constructor("QQ", ("x", "y"))

    unit_entry_case = _case(
        id = "ecp-unit-entry-gf2",
        kind = :unit_entry,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (a = x, b = y, c = one(R2)),
        column_order = (:a, :b, :c),
        monicity = _identity_monicity(:a, :x, x),
        witnesses = (),
        expected = (; current_status = :passes),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62",),
    )

    witness_unit_case = _case(
        id = "ecp-witness-unit-gf2",
        kind = :witness_unit,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (a = x, b = y, c = x + one(R2)),
        column_order = (:a, :b, :c),
        monicity = _identity_monicity(:a, :x, x),
        witnesses = ((
            kind = :ideal_membership,
            coefficients = (one(R2), zero(R2), one(R2)),
            require_nonunit_coefficients = false,
        ),),
        expected = (; current_status = :passes),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62",),
    )

    variable_change_case = _case(
        id = "ecp-variable-change-monic-gf2",
        kind = :variable_change,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (
            a = x + y^2,
            b = x * y + x + one(R2),
            c = x^2 + x * y + y + one(R2),
        ),
        column_order = (:a, :b, :c),
        monicity = (;
            selected_entry = :a,
            variable_name = :x,
            substitution = (; x = x + y),
            transformed_entry = x + y + y^2,
        ),
        witnesses = (),
        expected = (; current_status = :passes),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62",),
    )

    variable_change_permuted_case = _case(
        id = "ecp-variable-change-permuted-gf2",
        kind = :variable_change,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (
            a = x + y^2,
            b = x * y + x + one(R2),
            c = x^2 + x * y + y + one(R2),
        ),
        column_order = (:b, :a, :c),
        monicity = (;
            selected_entry = :b,
            variable_name = :y,
            substitution = (; x = x + y^2),
            transformed_entry = x * y + x + y^3 + y^2 + one(R2),
        ),
        witnesses = (),
        expected = (; current_status = :passes),
        source_refs = ("Issue 85 variable-change replay coverage",),
        consumer_issue_ids = ("#62", "#85"),
    )

    link_bezout_case = _case(
        id = "ecp-link-bezout-nonunit-witness-qq",
        kind = :link_bezout,
        stage_coverage = :supported,
        ring_constructor = qq_ring_constructor,
        ring = _ring_metadata("QQ[x, y]", RQ, ("x", "y"), (X, Y)),
        variable_order = (:x, :y),
        entries = (a = one(RQ) - X, b = X, c = X * Y),
        column_order = (:a, :b, :c),
        monicity = _identity_monicity(:b, :x, X),
        witnesses = ((
            kind = :link_bezout,
            coefficients = (one(RQ) + X, X, zero(RQ)),
            resultant = one(RQ),
            coverage = (; covers_unit_ideal = true),
            path = (:bezout, :coverage),
            require_nonunit_coefficients = true,
        ),),
        expected = (; current_status = :passes),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62", "#88"),
    )

    longer_embedded_block_case = _case(
        id = "ecp-longer-embedded-block-gf2",
        kind = :embedded_block,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (
            a = x + y^2,
            b = x * y + x + one(R2),
            c = x^2 + x * y + y + one(R2),
            d = x^2,
            e = x * y,
            f = y^2 + x,
        ),
        column_order = (:a, :b, :c, :d, :e, :f),
        monicity = _identity_monicity(:a, :x, x + y^2),
        witnesses = (),
        expected = (; current_status = :passes),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62",),
    )

    unsupported_case = _case(
        id = "ecp-unsupported-unimodular-gf2",
        kind = :unsupported,
        stage_coverage = :staged,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (a = zero(R2), b = x^2, c = x * y + one(R2)),
        column_order = (:a, :b, :c),
        monicity = _identity_monicity(:b, :x, x^2),
        witnesses = ((
            kind = :missing_link_witness,
            missing = (:exact_reduction_stage,),
        ),),
        expected = (; current_status = :staged_fail, message_substring = "unsupported exact unimodular column reduction"),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62", "#88"),
    )

    non_unimodular_case = _case(
        id = "ecp-non-unimodular-gf2",
        kind = :non_unimodular,
        stage_coverage = :unsupported,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (a = x, b = y, c = x * y),
        column_order = (:a, :b, :c),
        monicity = _identity_monicity(:a, :x, x),
        witnesses = (),
        expected = (; current_status = :rejects_non_unimodular, message_substring = "v must be a unimodular column"),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62",),
    )

    monic_first_entry_case = _case(
        id = "ecp-monic-first-entry-qq",
        kind = :monic_first_entry,
        stage_coverage = :supported,
        ring_constructor = qq_ring_constructor,
        ring = _ring_metadata("QQ[x, y]", RQ, ("x", "y"), (X, Y)),
        variable_order = (:x, :y),
        entries = (a = X^2 + Y + one(RQ), b = X, c = Y),
        column_order = (:a, :b, :c),
        monicity = _identity_monicity(:a, :x, X^2 + Y + one(RQ)),
        witnesses = (),
        expected = (; current_status = :passes),
        source_refs = ("Suslin exact unimodular reduction tests",),
        consumer_issue_ids = ("#62", "#87"),
    )

    bad_witness = _negative_control(
        "ecp-corrupt-witness-control",
        "ecp-witness-unit-gf2",
        "mutated witness coefficients",
        merge(
            witness_unit_case,
            (;
                witnesses = (merge(only(witness_unit_case.witnesses), (;
                    coefficients = (zero(R2), zero(R2), one(R2)),
                )),),
            ),
        ),
    )

    bad_monicity = _negative_control(
        "ecp-corrupt-monicity-control",
        "ecp-variable-change-monic-gf2",
        "mutated monicity transformed entry",
        merge(
            variable_change_case,
            (;
                monicity = merge(variable_change_case.monicity, (;
                    transformed_entry = y,
                )),
            ),
        ),
    )

    return (;
        cases = [
            unit_entry_case,
            witness_unit_case,
            variable_change_case,
            variable_change_permuted_case,
            link_bezout_case,
            longer_embedded_block_case,
            unsupported_case,
            non_unimodular_case,
            monic_first_entry_case,
        ],
        negative_controls = [
            bad_witness,
            bad_monicity,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
