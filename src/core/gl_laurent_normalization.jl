function _laurent_monomial_metadata(value)
    try
        exponent_vectors = collect(exponents(value))
        coeffs = collect(coefficients(value))
        length(exponent_vectors) == 1 || return nothing
        length(coeffs) == 1 || return nothing
        return (;
            monomial_exponents = Tuple(Int.(collect(exponent_vectors[1]))),
            monomial_coefficient = coeffs[1],
        )
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        return nothing
    end
end

function _laurent_monomial_term(R, coeff, exponent_vector)
    term = R(coeff)
    for (j, exponent) in enumerate(exponent_vector)
        exponent == 0 && continue
        term *= gen(R, j)^exponent
    end
    return term
end

function _is_laurent_other_unit_via_nilpotent_perturbation(determinant, R)::Bool
    exponent_vectors = collect(exponents(determinant))
    coeffs = collect(coefficients(determinant))
    length(exponent_vectors) == length(coeffs) || return false
    length(coeffs) > 1 || return false

    for (coeff, raw_exponents) in zip(coeffs, exponent_vectors)
        coeff_is_unit = try
            is_unit(coeff)
        catch err
            err isa MethodError || rethrow()
            false
        end
        coeff_is_unit || continue

        exponent_vector = Int.(collect(raw_exponents))
        term = _laurent_monomial_term(R, coeff, exponent_vector)
        is_nilpotent(inv(term) * determinant - one(R)) && return true
    end

    return false
end

function _is_supported_laurent_monomial_unit(monomial)::Bool
    try
        return is_unit(monomial.monomial_coefficient)
    catch err
        err isa MethodError || rethrow()
        return false
    end
end

function classify_laurent_determinant(A)
    nrows(A) == ncols(A) || throw(ArgumentError("Laurent GL_n determinant classification requires a square matrix"))
    R = _require_laurent_polynomial_ring(base_ring(A); label="A base ring")
    determinant = det(A)
    determinant == one(R) && return (;
        determinant,
        classification = :one,
        monomial_exponents = ntuple(_ -> 0, ngens(R)),
        monomial_coefficient = one(R),
    )

    negative_one = -one(R)
    if negative_one != one(R) && determinant == negative_one
        return (;
            determinant,
            classification = :permutation_sign_unit,
            monomial_exponents = ntuple(_ -> 0, ngens(R)),
            monomial_coefficient = negative_one,
        )
    end

    monomial = _laurent_monomial_metadata(determinant)
    if monomial !== nothing
        if _is_supported_laurent_monomial_unit(monomial)
            return (;
                determinant,
                classification = :laurent_monomial_unit,
                monomial_exponents = monomial.monomial_exponents,
                monomial_coefficient = monomial.monomial_coefficient,
            )
        end

        return (;
            determinant,
            classification = :non_unit,
            monomial_exponents = nothing,
            monomial_coefficient = nothing,
        )
    end

    _is_laurent_other_unit_via_nilpotent_perturbation(determinant, R) && return (;
        determinant,
        classification = :other_unit,
        monomial_exponents = nothing,
        monomial_coefficient = nothing,
    )

    determinant_is_unit = try
        is_unit(determinant)
    catch err
        err isa MethodError || rethrow()
        false
    end

    determinant_is_unit || return (;
        determinant,
        classification = :non_unit,
        monomial_exponents = nothing,
        monomial_coefficient = nothing,
    )

    return (;
        determinant,
        classification = :other_unit,
        monomial_exponents = nothing,
        monomial_coefficient = nothing,
    )
end

function _identity_correction(R, n::Int, determinant)
    identity = identity_matrix(R, n)
    return (;
        kind = :identity,
        side = :left,
        factor = identity,
        inverse_factor = identity,
        determinant,
    )
end

function _expected_laurent_correction_kind(classification::Symbol)
    return classification == :one ? :identity : :left_diagonal_determinant_correction
end

function _left_diagonal_determinant_correction(R, n::Int, determinant)
    factor = identity_matrix(R, n)
    factor[1, 1] = determinant
    inverse_factor = identity_matrix(R, n)
    inverse_factor[1, 1] = inv(determinant)
    return (;
        kind = :left_diagonal_determinant_correction,
        side = :left,
        factor,
        inverse_factor,
        determinant,
    )
end

function _throw_unsupported_laurent_gl_determinant(classification)
    if classification == :non_unit
        throw(ArgumentError("unsupported Laurent GL_n determinant: determinant is non-unit, so the input is outside the staged SL_n factorization path"))
    elseif classification == :other_unit
        throw(ArgumentError("unsupported Laurent GL_n determinant: non-monomial units are outside the staged SL_n factorization path"))
    end

    return nothing
end

function normalize_laurent_gl_matrix(A)
    nrows(A) == ncols(A) || throw(ArgumentError("Laurent GL_n normalization requires a square matrix"))
    R = _require_laurent_polynomial_ring(base_ring(A); label="A base ring")
    n = nrows(A)
    determinant_profile = classify_laurent_determinant(A)
    classification = determinant_profile.classification

    _throw_unsupported_laurent_gl_determinant(classification)

    correction = classification == :one ?
        _identity_correction(R, n, determinant_profile.determinant) :
        _left_diagonal_determinant_correction(R, n, determinant_profile.determinant)
    normalized_matrix = correction.inverse_factor * A
    normalization = (;
        input_size = (n, n),
        ring = R,
        determinant = determinant_profile.determinant,
        determinant_classification = classification,
        determinant_profile,
        normalized_matrix,
        correction,
    )

    verify_laurent_gl_normalization(A, normalization) || throw(ArgumentError("Laurent GL_n normalization failed exact reconstruction verification"))
    return normalization
end

function verify_laurent_gl_normalization(A, normalization)::Bool
    try
        nrows(A) == ncols(A) || return false
        n = nrows(A)
        R = base_ring(A)
        _require_laurent_polynomial_ring(R; label="A base ring")
        normalization.input_size == (n, n) || return false
        (normalization.ring === R || normalization.ring == R) || return false
        determinant_profile = classify_laurent_determinant(A)
        normalization.determinant == determinant_profile.determinant || return false
        normalization.determinant_classification == determinant_profile.classification || return false
        nrows(normalization.normalized_matrix) == n || return false
        ncols(normalization.normalized_matrix) == n || return false
        base_ring(normalization.normalized_matrix) == R || return false
        correction = normalization.correction
        correction.side == :left || return false
        correction.kind == _expected_laurent_correction_kind(normalization.determinant_classification) || return false
        correction.determinant == normalization.determinant || return false
        nrows(correction.factor) == n || return false
        ncols(correction.factor) == n || return false
        nrows(correction.inverse_factor) == n || return false
        ncols(correction.inverse_factor) == n || return false
        base_ring(correction.factor) == R || return false
        base_ring(correction.inverse_factor) == R || return false
        identity = identity_matrix(R, n)
        correction.factor * correction.inverse_factor == identity || return false
        correction.inverse_factor * correction.factor == identity || return false
        correction.factor * normalization.normalized_matrix == A || return false
        det(normalization.normalized_matrix) == one(R) || return false
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
