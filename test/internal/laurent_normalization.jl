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
    vector = [column[i, 1] for i in 1:nrows(column)]
    vector_normalization = normalize_laurent_object(vector)
    @test vector_normalization.metadata.kind == :vector
    @test vector_normalization.metadata.shape == (2,)
    @test vector_normalization.metadata.column_shifts == ((2, 0),)
    @test _all_polynomial_vector_entries_have_nonnegative_exponents(vector_normalization.normalized_object)
    lifted_vector = lift_laurent_normalization(vector_normalization)
    @test eltype(lifted_vector) == typeof(first(vector))
    @test lifted_vector == vector

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

    laurent_normalized_column = only(column_normalization.metadata.shift_monomials) * column
    laurent_ring_metadata = merge(column_normalization.metadata, (; polynomial_ring = R))
    laurent_ring_tampered = (;
        normalized_object = laurent_normalized_column,
        metadata = laurent_ring_metadata,
    )
    @test_throws ArgumentError lift_laurent_normalization(laurent_ring_tampered)
    @test !verify_laurent_normalization(column, laurent_ring_tampered)

    invalid_polynomial_ring_metadata = merge(column_normalization.metadata, (; polynomial_ring = 42))
    invalid_polynomial_ring_tampered = (;
        normalized_object = column_normalization.normalized_object,
        metadata = invalid_polynomial_ring_metadata,
    )
    @test_throws ArgumentError lift_laurent_normalization(invalid_polynomial_ring_tampered)

    @test_throws ArgumentError normalize_laurent_object(one(R))

    determinant_on_column_metadata = merge(
        column_normalization.metadata,
        (; determinant_shift_exponents = (2, 0)),
    )
    determinant_on_column_tampered = (;
        normalized_object = column_normalization.normalized_object,
        metadata = determinant_on_column_metadata,
    )
    @test_throws ArgumentError lift_laurent_normalization(determinant_on_column_tampered)

    unsupported_kind_metadata = merge(vector_normalization.metadata, (; kind = :scalar))
    unsupported_kind_tampered = (;
        normalized_object = vector_normalization.normalized_object,
        metadata = unsupported_kind_metadata,
    )
    @test_throws ArgumentError lift_laurent_normalization(unsupported_kind_tampered)
end
