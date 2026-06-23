module ParkWoodburnPolynomialFixtureCatalog

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
        role,
        route,
        status,
        provenance,
        ring_constructor,
        ring,
        matrix,
        determinant_expectation,
        source_refs,
        consumer_issue_ids)
    return (;
        id = id,
        role = role,
        route = route,
        status = status,
        provenance = provenance,
        ring_constructor = ring_constructor,
        ring = ring,
        matrix = matrix,
        determinant_expectation = determinant_expectation,
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

function _product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function catalog()
    RQ, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    R2, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    RR, (PX, Pr, Pg) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])

    qq_ring_constructor = _ordinary_ring_constructor("QQ", ("X",))
    gf2_ring_constructor = _ordinary_ring_constructor("GF(2)", ("x", "y"))
    quillen_ring_constructor = _ordinary_ring_constructor("QQ", ("X", "r", "g"))

    qq_ring = _ring_metadata("QQ[X]", RQ, ("X",), (X,))
    gf2_ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y))
    quillen_ring = _ring_metadata("QQ[X, r, g]", RR, ("X", "r", "g"), (PX, Pr, Pg))

    fast_local_matrix = matrix(RQ, [
        one(RQ) X + one(RQ) zero(RQ);
        X one(RQ) + X + X^2 zero(RQ);
        zero(RQ) zero(RQ) one(RQ)
    ])
    fast_local_case = _case(
        id = "pw-poly-univariate-sl3-fast-local-qq",
        role = :univariate_sl3_fast_local,
        route = :fast_local_sl3,
        status = :supported,
        provenance = (;
            source = "Park-Woodburn issue 109 univariate SL3 fast-local shell",
        ),
        ring_constructor = qq_ring_constructor,
        ring = qq_ring,
        matrix = fast_local_matrix,
        determinant_expectation = :one,
        source_refs = ("Park-Woodburn issue 109 fast-local SL_3 witness",),
        consumer_issue_ids = ("#64", "#109", "#110"),
    )

    block_a = matrix(RQ, [
        one(RQ) + X one(RQ) zero(RQ);
        X one(RQ) zero(RQ);
        zero(RQ) zero(RQ) one(RQ)
    ])
    block_b = matrix(RQ, [
        one(RQ) one(RQ) + X zero(RQ);
        X one(RQ) + X + X^2 zero(RQ);
        zero(RQ) zero(RQ) one(RQ)
    ])
    disjoint_blocks_matrix = matrix(RQ, [
        one(RQ) + X one(RQ) zero(RQ) zero(RQ) zero(RQ) zero(RQ);
        X one(RQ) zero(RQ) zero(RQ) zero(RQ) zero(RQ);
        zero(RQ) zero(RQ) one(RQ) zero(RQ) zero(RQ) zero(RQ);
        zero(RQ) zero(RQ) zero(RQ) one(RQ) one(RQ) + X zero(RQ);
        zero(RQ) zero(RQ) zero(RQ) X one(RQ) + X + X^2 zero(RQ);
        zero(RQ) zero(RQ) zero(RQ) zero(RQ) zero(RQ) one(RQ)
    ])
    disjoint_blocks_case = _case(
        id = "pw-poly-univariate-sln-disjoint-blocks-qq",
        role = :univariate_sln_disjoint_blocks,
        route = :disjoint_local_blocks,
        status = :supported,
        provenance = (;
            source = "Park-Woodburn issue 109 disjoint local SL3 blocks",
        ),
        ring_constructor = qq_ring_constructor,
        ring = qq_ring,
        matrix = disjoint_blocks_matrix,
        determinant_expectation = :one,
        source_refs = ("test/expert/sln_to_sl3_reduction.jl local block product",),
        consumer_issue_ids = ("#64", "#109", "#110"),
    )

    recursive_factors = (
        elementary_matrix(4, 1, 2, x + y, R2),
        elementary_matrix(4, 2, 3, x * y + x, R2),
        elementary_matrix(4, 3, 4, y + one(R2), R2),
    )
    recursive_matrix = _product(recursive_factors, R2, 4)
    recursive_case = _case(
        id = "pw-poly-recursive-column-peel-gf2",
        role = :recursive_column_peel,
        route = :recursive_column_peel,
        status = :staged,
        provenance = (;
            source = "Park-Woodburn issue 109 recursive column peel staging example",
        ),
        ring_constructor = gf2_ring_constructor,
        ring = gf2_ring,
        matrix = recursive_matrix,
        determinant_expectation = :one,
        source_refs = ("Park-Woodburn issue 109 recursive peel staging matrix",),
        consumer_issue_ids = ("#64", "#109", "#113"),
    )

    quillen_matrix = _product(
        (
            elementary_matrix(3, 1, 2, PX + Pr + Pg, RR),
            elementary_matrix(3, 2, 3, Pr * Pg + PX, RR),
        ),
        RR,
        3,
    )
    quillen_case = _case(
        id = "quillen-patched-substitution-witness-qq",
        role = :multivariate_quillen,
        route = :quillen_patched_substitution,
        status = :blocked,
        provenance = (;
            source = "Issue 99 quillen patched substitution witness",
            quillen_fixture_id = "quillen-patched-substitution-witness-qq",
        ),
        ring_constructor = quillen_ring_constructor,
        ring = quillen_ring,
        matrix = quillen_matrix,
        determinant_expectation = :one,
        source_refs = ("Issue 99 quillen patched substitution witness shape",),
        consumer_issue_ids = ("#64", "#99", "#105", "#109", "#115"),
    )

    det_not_one_matrix = matrix(RQ, [
        one(RQ) + X one(RQ) zero(RQ);
        zero(RQ) one(RQ) zero(RQ);
        zero(RQ) zero(RQ) one(RQ)
    ])
    det_not_one_control = _negative_control(
        "pw-poly-det-not-one-control",
        "pw-poly-univariate-sl3-fast-local-qq",
        "matrix determinant is not one despite determinant-one metadata",
        merge(fast_local_case, (;
            matrix = det_not_one_matrix,
        )),
    )

    outside_witness_control = _negative_control(
        "pw-poly-det-one-outside-witness-control",
        "pw-poly-univariate-sl3-fast-local-qq",
        "determinant-one matrix claims an unsupported route",
        merge(fast_local_case, (;
            route = :unsupported_route_family,
            status = :supported,
            provenance = (;
                source = "negative control for unsupported route metadata",
            ),
        )),
    )

    wrong_route_control = _negative_control(
        "pw-poly-wrong-route-control",
        "pw-poly-univariate-sl3-fast-local-qq",
        "fast-local matrix is mislabeled as another currently supported route",
        merge(fast_local_case, (;
            route = :disjoint_local_blocks,
            status = :supported,
        )),
    )

    return (;
        cases = [
            fast_local_case,
            disjoint_blocks_case,
            recursive_case,
            quillen_case,
        ],
        negative_controls = [
            det_not_one_control,
            outside_witness_control,
            wrong_route_control,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
