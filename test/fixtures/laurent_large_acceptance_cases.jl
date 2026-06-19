module LaurentLargeAcceptanceCases

using Oscar
using Suslin

include("toricbuilder_factor_toric_block_3.jl")

const ACCEPTANCE_40_SIZE = 40
const ACCEPTANCE_LARGER_SIZE = 48

function _block_locations(n::Int)
    return [Int[first_idx, first_idx + 1, first_idx + 2] for first_idx in 1:3:(n - 2)]
end

function _local_laurent_block(R, x, y, block_index::Int)
    q = isodd(block_index) ? x * y : x + y^-1
    r = one(R)
    p = one(R) + q * r
    s = one(R)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _large_block_local_case(R, x, y, n::Int, id::String, description::String)
    locations = _block_locations(n)
    A = identity_matrix(R, n)
    for (block_index, indices) in enumerate(locations)
        A *= block_embedding(_local_laurent_block(R, x, y, block_index), n, indices)
    end

    return (;
        id,
        kind = :large_laurent_factorization,
        ring = (;
            description = "GF(2)[x^+/-1, y^+/-1]",
            object = R,
            generators = (x, y),
        ),
        size = (n, n),
        matrix = A,
        block_locations = locations,
        expected_path = :elementary_factorization,
        provenance = (;
            source = :synthetic_supported_block_local,
            issue = "#17",
            description,
            construction = "product of disjoint embedded Laurent local SL3 blocks",
        ),
        negative_control = (;
            kind = :replace_first_factor_with_identity,
            description = "Replacing the first returned factor by identity must break exact reconstruction.",
        ),
    )
end

function _toricbuilder_pinv_case()
    fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    pinv = only(filter(entry -> entry.toricbuilder_role == "Pinv", fixture.cases))
    return (;
        id = "toricbuilder-factor-toric-block-3-pinv-normalized-contract",
        kind = :toricbuilder_normalized_contract,
        ring = pinv.ring,
        size = pinv.size,
        matrix = pinv.matrix,
        source_matrix = pinv.source_matrix,
        expected_path = :normalized_contract,
        provenance = (;
            source = :toricbuilder_contract_fixture,
            issue = "#19",
            fixture_id = pinv.name,
            toricbuilder_role = pinv.toricbuilder_role,
            toricbuilder_commit = pinv.provenance.toricbuilder_commit,
            generation_command = pinv.provenance.generation_command,
        ),
        negative_control = (;
            kind = :corrupt_inverse_relation_entry,
            row = 1,
            col = 1,
            description = "Toggling one Pinv entry must break the exact ToricBuilder inverse relation.",
        ),
    )
end

function acceptance_catalog()
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    return (;
        cases = [
            _toricbuilder_pinv_case(),
            _large_block_local_case(
                R,
                x,
                y,
                ACCEPTANCE_40_SIZE,
                "laurent-block-local-40x40",
                "40x40 determinant-one Laurent matrix with disjoint supported local SL3 blocks",
            ),
            _large_block_local_case(
                R,
                x,
                y,
                ACCEPTANCE_LARGER_SIZE,
                "laurent-block-local-48x48",
                "48x48 determinant-one Laurent matrix with disjoint supported local SL3 blocks",
            ),
        ],
    )
end

end
