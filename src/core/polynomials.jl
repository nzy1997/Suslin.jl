# Reserved split point for polynomial helpers that the later algorithm layers build on.

function _laurent_normalization_variable_names(R)
    return Tuple(string.(symbols(R)))
end

function _ordinary_polynomial_ring_for_laurent(R)
    P, _ = suslin_polynomial_ring(coefficient_ring(R), collect(_laurent_normalization_variable_names(R)))
    return P
end

function _laurent_monomial_from_exponents(R, exponent_vector)
    length(exponent_vector) == ngens(R) || throw(ArgumentError("shift exponent vector length must match the Laurent ring generators"))
    term = one(R)
    for (j, exponent) in enumerate(exponent_vector)
        exponent == 0 && continue
        term *= gen(R, j)^exponent
    end
    return term
end

function _column_shift_for_laurent_entries(entries, R)
    min_exponents = zeros(Int, ngens(R))
    for value in entries
        _require_laurent_element(value, R; label="normalization entry")
        for raw_exponents in collect(exponents(value))
            exponent_vector = Int.(collect(raw_exponents))
            length(exponent_vector) == ngens(R) || throw(ArgumentError("Laurent exponent vector length must match the parent ring"))
            for j in 1:ngens(R)
                min_exponents[j] = min(min_exponents[j], exponent_vector[j])
            end
        end
    end
    return ntuple(j -> max(0, -min_exponents[j]), ngens(R))
end

function _polynomial_term_from_exponents(P, coeff, exponent_vector)
    term = P(coeff)
    for (j, exponent) in enumerate(exponent_vector)
        exponent < 0 && throw(ArgumentError("normalized Laurent term still has a negative exponent"))
        exponent == 0 && continue
        term *= gen(P, j)^exponent
    end
    return term
end

function _laurent_entry_to_polynomial(value, R, P, shift)
    _require_laurent_element(value, R; label="normalization entry")
    length(shift) == ngens(R) || throw(ArgumentError("shift exponent vector length must match the Laurent ring generators"))

    result = zero(P)
    for (coeff, raw_exponents) in zip(collect(coefficients(value)), collect(exponents(value)))
        exponent_vector = Int.(collect(raw_exponents))
        shifted_exponents = ntuple(j -> exponent_vector[j] + shift[j], ngens(R))
        result += _polynomial_term_from_exponents(P, coeff, shifted_exponents)
    end
    return result
end

function _sum_laurent_shift_exponents(column_shifts, R)
    return ntuple(j -> sum(shift[j] for shift in column_shifts), ngens(R))
end

function _laurent_normalization_metadata(kind::Symbol, shape, R, P, column_shifts)
    shift_monomials = map(shift -> _laurent_monomial_from_exponents(R, shift), column_shifts)
    inverse_shift_monomials = map(
        shift -> _laurent_monomial_from_exponents(R, ntuple(j -> -shift[j], ngens(R))),
        column_shifts,
    )
    determinant_shift_exponents = kind == :matrix && shape[1] == shape[2] ?
        _sum_laurent_shift_exponents(column_shifts, R) :
        nothing
    return (;
        kind,
        shape,
        laurent_ring = R,
        polynomial_ring = P,
        variable_names = _laurent_normalization_variable_names(R),
        column_shifts,
        shift_monomials,
        inverse_shift_monomials,
        determinant_shift_exponents,
    )
end

function _has_oscar_matrix_interface(obj)::Bool
    try
        nrows(obj)
        ncols(obj)
        base_ring(obj)
        return true
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        return false
    end
end

function _normalize_laurent_matrix(A)
    R = _require_laurent_polynomial_ring(base_ring(A); label="input base ring")
    rows, cols = nrows(A), ncols(A)
    P = _ordinary_polynomial_ring_for_laurent(R)
    column_shifts = ntuple(
        j -> _column_shift_for_laurent_entries((A[i, j] for i in 1:rows), R),
        cols,
    )

    normalized = zero_matrix(P, rows, cols)
    for j in 1:cols
        for i in 1:rows
            normalized[i, j] = _laurent_entry_to_polynomial(A[i, j], R, P, column_shifts[j])
        end
    end

    metadata = _laurent_normalization_metadata(:matrix, (rows, cols), R, P, column_shifts)
    return (; normalized_object = normalized, metadata)
end

function _normalize_laurent_vector(values::AbstractVector)
    R = _require_same_laurent_parent(values; label="vector entries")
    P = _ordinary_polynomial_ring_for_laurent(R)
    shift = _column_shift_for_laurent_entries(values, R)
    normalized = [_laurent_entry_to_polynomial(value, R, P, shift) for value in values]
    metadata = _laurent_normalization_metadata(:vector, (length(values),), R, P, (shift,))
    return (; normalized_object = normalized, metadata)
end

function normalize_laurent_object(obj)
    _has_oscar_matrix_interface(obj) && return _normalize_laurent_matrix(obj)
    obj isa AbstractVector && return _normalize_laurent_vector(obj)
    throw(ArgumentError("input must be an Oscar matrix over a Laurent polynomial ring or a vector of Laurent elements"))
end

function _require_laurent_normalization_metadata(metadata)
    required = (
        :kind,
        :shape,
        :laurent_ring,
        :polynomial_ring,
        :variable_names,
        :column_shifts,
        :shift_monomials,
        :inverse_shift_monomials,
        :determinant_shift_exponents,
    )
    for field in required
        hasproperty(metadata, field) || throw(ArgumentError("normalization metadata missing field $(field)"))
    end

    R = _require_laurent_polynomial_ring(metadata.laurent_ring; label="metadata Laurent ring")
    metadata.variable_names == _laurent_normalization_variable_names(R) || throw(ArgumentError("normalization metadata variable names do not match the Laurent ring"))
    metadata.column_shifts isa Tuple || throw(ArgumentError("normalization metadata column shifts must be a tuple"))

    for shift in metadata.column_shifts
        length(shift) == ngens(R) || throw(ArgumentError("shift exponent vector length must match the Laurent ring generators"))
        all(exponent -> exponent >= 0, shift) || throw(ArgumentError("shift exponents must be nonnegative"))
    end

    expected_shift_monomials = map(shift -> _laurent_monomial_from_exponents(R, shift), metadata.column_shifts)
    expected_inverse_shift_monomials = map(
        shift -> _laurent_monomial_from_exponents(R, ntuple(j -> -shift[j], ngens(R))),
        metadata.column_shifts,
    )
    metadata.shift_monomials == expected_shift_monomials || throw(ArgumentError("normalization shift monomials do not match shift exponents"))
    metadata.inverse_shift_monomials == expected_inverse_shift_monomials || throw(ArgumentError("normalization inverse shift monomials do not match shift exponents"))

    if metadata.kind == :matrix && metadata.shape[1] == metadata.shape[2]
        metadata.determinant_shift_exponents == _sum_laurent_shift_exponents(metadata.column_shifts, R) ||
            throw(ArgumentError("normalization determinant shift metadata does not match column shifts"))
    elseif metadata.determinant_shift_exponents !== nothing
        throw(ArgumentError("normalization determinant shift metadata is only defined for square matrices"))
    end

    return metadata
end

function _lift_laurent_matrix(polynomial_object, metadata)
    R = metadata.laurent_ring
    P = metadata.polynomial_ring
    rows, cols = metadata.shape
    nrows(polynomial_object) == rows || throw(ArgumentError("normalized matrix row count does not match metadata"))
    ncols(polynomial_object) == cols || throw(ArgumentError("normalized matrix column count does not match metadata"))
    base_ring(polynomial_object) == P || throw(ArgumentError("normalized matrix base ring does not match metadata"))
    length(metadata.column_shifts) == cols || throw(ArgumentError("matrix metadata must record one shift per column"))

    lifted = zero_matrix(R, rows, cols)
    for j in 1:cols
        inverse_shift = metadata.inverse_shift_monomials[j]
        for i in 1:rows
            lifted[i, j] = inverse_shift * R(polynomial_object[i, j])
        end
    end
    return lifted
end

function _lift_laurent_vector(polynomial_object, metadata)
    R = metadata.laurent_ring
    P = metadata.polynomial_ring
    length(polynomial_object) == metadata.shape[1] || throw(ArgumentError("normalized vector length does not match metadata"))
    length(metadata.column_shifts) == 1 || throw(ArgumentError("vector metadata must record exactly one shift"))
    inverse_shift = only(metadata.inverse_shift_monomials)

    lifted = Vector{Any}(undef, length(polynomial_object))
    for (i, value) in enumerate(polynomial_object)
        parent(value) == P || throw(ArgumentError("normalized vector entry parent does not match metadata"))
        lifted[i] = inverse_shift * R(value)
    end
    return lifted
end

function lift_laurent_normalization(polynomial_object, metadata)
    metadata = _require_laurent_normalization_metadata(metadata)
    if metadata.kind == :matrix
        return _lift_laurent_matrix(polynomial_object, metadata)
    elseif metadata.kind == :vector
        return _lift_laurent_vector(polynomial_object, metadata)
    end
    throw(ArgumentError("unsupported Laurent normalization kind $(metadata.kind)"))
end

function lift_laurent_normalization(normalization)
    hasproperty(normalization, :normalized_object) || throw(ArgumentError("normalization missing normalized_object"))
    hasproperty(normalization, :metadata) || throw(ArgumentError("normalization missing metadata"))
    return lift_laurent_normalization(normalization.normalized_object, normalization.metadata)
end

function verify_laurent_normalization(original, normalization)::Bool
    try
        return lift_laurent_normalization(normalization) == original
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
