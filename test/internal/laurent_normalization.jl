using Test
using Suslin
using Oscar

const LAURENT_NORMALIZATION_FIXTURES = joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl")

function _fixture_by_id(catalog, id)
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _all_polynomial_matrix_entries_have_nonnegative_exponents(A)
    for value in A
        for exponent_vector in collect(exponents(value))
            all(exponent -> exponent >= 0, exponent_vector) || return false
        end
    end
    return true
end

function _all_polynomial_vector_entries_have_nonnegative_exponents(values)
    for value in values
        for exponent_vector in collect(exponents(value))
            all(exponent -> exponent >= 0, exponent_vector) || return false
        end
    end
    return true
end

@testset "Laurent normalization and lift-back helpers" begin
    include(LAURENT_NORMALIZATION_FIXTURES)
    catalog = LaurentFixtureCatalog.catalog()

    column_fixture = _fixture_by_id(catalog, "laurent-negative-exponent-normalization")
    column = column_fixture.inputs.vector
    column_normalization = normalize_laurent_object(column)

    @test column_normalization.metadata.kind == :matrix
    @test column_normalization.metadata.shape == (2, 1)
    @test column_normalization.metadata.column_shifts == ((2, 0),)
    @test column_normalization.metadata.determinant_shift_exponents === nothing
    @test base_ring(column_normalization.normalized_object) == column_normalization.metadata.polynomial_ring
    @test _all_polynomial_matrix_entries_have_nonnegative_exponents(column_normalization.normalized_object)
    @test lift_laurent_normalization(column_normalization) == column
    @test lift_laurent_normalization(column_normalization.normalized_object, column_normalization.metadata) == column
    @test verify_laurent_normalization(column, column_normalization)

    matrix_fixture = _fixture_by_id(catalog, "toricbuilder-factor-toric-block-3-pinv")
    A = matrix_fixture.inputs.matrix
    matrix_normalization = normalize_laurent_object(A)

    @test matrix_normalization.metadata.kind == :matrix
    @test matrix_normalization.metadata.shape == (8, 8)
    @test length(matrix_normalization.metadata.column_shifts) == 8
    @test matrix_normalization.metadata.determinant_shift_exponents ==
        ntuple(i -> sum(shift[i] for shift in matrix_normalization.metadata.column_shifts), 2)
    @test base_ring(matrix_normalization.normalized_object) == matrix_normalization.metadata.polynomial_ring
    @test _all_polynomial_matrix_entries_have_nonnegative_exponents(matrix_normalization.normalized_object)
    @test lift_laurent_normalization(matrix_normalization) == A
    @test verify_laurent_normalization(A, matrix_normalization)

    R = column_fixture.ring.object
    x, y = column_fixture.ring.generators
    vector_normalization = normalize_laurent_object([x^-1, x^-2 * y])
    @test vector_normalization.metadata.kind == :vector
    @test vector_normalization.metadata.shape == (2,)
    @test vector_normalization.metadata.column_shifts == ((2, 0),)
    @test _all_polynomial_vector_entries_have_nonnegative_exponents(vector_normalization.normalized_object)
    @test lift_laurent_normalization(vector_normalization) == [x^-1, x^-2 * y]

    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    @test_throws ArgumentError normalize_laurent_object([x, u])

    tampered_shift = ((1, 0),)
    tampered_metadata = merge(
        column_normalization.metadata,
        (;
            column_shifts = tampered_shift,
            shift_monomials = (x,),
            inverse_shift_monomials = (x^-1,),
        ),
    )
    tampered = (;
        normalized_object = column_normalization.normalized_object,
        metadata = tampered_metadata,
    )
    @test lift_laurent_normalization(tampered) != column
    @test !verify_laurent_normalization(column, tampered)
end
