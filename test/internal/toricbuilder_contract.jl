using Test
using Suslin
using Oscar

const TORICBUILDER_FIXTURE_PATH = joinpath(@__DIR__, "..", "fixtures", "toricbuilder_factor_toric_block_3.jl")

function determinant_classification(A)
    R = base_ring(A)
    d = det(A)
    d == one(R) && return "one"
    is_unit(d) && return "laurent_monomial_unit"
    return "non-unit"
end

function assert_required_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return true
end

function assert_toricbuilder_contract_entry_valid(entry)
    for field in (
        :name,
        :toricbuilder_role,
        :matrix,
        :source_matrix,
        :relation_description,
        :ring,
        :size,
        :determinant_classification,
        :expected_suslin_status,
        :expected_output,
        :provenance,
    )
        assert_required_field(entry, field)
    end

    nrows(entry.matrix) == entry.size[1] || throw(ArgumentError("fixture $(entry.name) row count does not match metadata"))
    ncols(entry.matrix) == entry.size[2] || throw(ArgumentError("fixture $(entry.name) column count does not match metadata"))
    entry.size[1] == entry.size[2] || throw(ArgumentError("fixture $(entry.name) must be square"))
    base_ring(entry.matrix) == base_ring(entry.source_matrix) || throw(ArgumentError("fixture $(entry.name) relation matrices use different rings"))
    entry.ring == "GF(2)[x^+/-1, y^+/-1]" || throw(ArgumentError("fixture $(entry.name) has unexpected ring metadata"))

    actual_class = determinant_classification(entry.matrix)
    actual_class == entry.determinant_classification || throw(ArgumentError("fixture $(entry.name) determinant classification $(entry.determinant_classification) does not match $(actual_class)"))

    R = base_ring(entry.matrix)
    expected_identity = identity_matrix(R, entry.size[1])
    entry.source_matrix * entry.matrix == expected_identity || throw(ArgumentError("fixture $(entry.name) exact ToricBuilder relation failed"))
    return true
end

@testset "toricbuilder contract fixture" begin
    @test isfile(TORICBUILDER_FIXTURE_PATH)

    if isfile(TORICBUILDER_FIXTURE_PATH)
        include(TORICBUILDER_FIXTURE_PATH)
        fixture = ToricBuilderFactorToricBlock3Fixture.fixture()

        @test fixture.ring.description == "GF(2)[x^+/-1, y^+/-1]"
        @test fixture.block_size == (2, 2)
        @test fixture.provenance.toricbuilder_commit == "fa7f82252d42fdc0b2726bc48af24ac4c70a8d73"

        roles = Set(entry.toricbuilder_role for entry in fixture.cases)
        @test "Qinv" in roles
        @test "Pinv" in roles

        for entry in fixture.cases
            @test assert_toricbuilder_contract_entry_valid(entry)
            @test entry.expected_output == :verified_transformation_certificate

            factors = elementary_factorization(entry.matrix)
            @test !isempty(factors)
            @test verify_factorization(entry.matrix, factors)
        end

        qinv = only(filter(entry -> entry.toricbuilder_role == "Qinv", fixture.cases))
        bad_class = merge(qinv, (; determinant_classification = "non-unit"))
        @test_throws ArgumentError assert_toricbuilder_contract_entry_valid(bad_class)

        corrupted_source = copy(qinv.source_matrix)
        corrupted_source[1, 1] += one(base_ring(corrupted_source))
        bad_relation = merge(qinv, (; source_matrix = corrupted_source))
        @test_throws ArgumentError assert_toricbuilder_contract_entry_valid(bad_relation)
    end
end
