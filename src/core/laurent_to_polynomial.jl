struct LaurentNoetherCertificate
    original_column
    ring
    selected_entry_index::Int
    selected_generator_index::Int
    selected_generator
    other_generator_index::Int
    noether_power::Int
    forward_substitution
    inverse_substitution
    transformed_column
    replayed_selected_entry
    leading_coefficient
    trailing_coefficient
    leading_coefficient_is_unit::Bool
    trailing_coefficient_is_unit::Bool
    validation_status::Symbol
end

function _laurent_noether_substitution_values(substitution_map)
    try
        return tuple((entry.value for entry in substitution_map)...)
    catch err
        err isa InterruptException && rethrow()
        return ()
    end
end

function _laurent_noether_map(variables, values)
    return tuple(((; variable = variables[index], value = values[index]) for index in eachindex(variables))...)
end

function _laurent_noether_maps_are_inverse(R, forward_substitution, inverse_substitution)::Bool
    try
        variables = tuple(gens(R)...)
        forward_values = _laurent_noether_substitution_values(forward_substitution)
        inverse_values = _laurent_noether_substitution_values(inverse_substitution)
        length(forward_values) == length(variables) || return false
        length(inverse_values) == length(variables) || return false
        all(entry -> hasproperty(entry, :variable) && hasproperty(entry, :value), forward_substitution) || return false
        all(entry -> hasproperty(entry, :variable) && hasproperty(entry, :value), inverse_substitution) || return false
        tuple((entry.variable for entry in forward_substitution)...) == variables || return false
        tuple((entry.variable for entry in inverse_substitution)...) == variables || return false
        forward_then_inverse = tuple((
            _coerce_into_ring(R, evaluate(evaluate(variable, collect(forward_values)), collect(inverse_values)), "forward-inverse Laurent substitution")
            for variable in variables
        )...)
        inverse_then_forward = tuple((
            _coerce_into_ring(R, evaluate(evaluate(variable, collect(inverse_values)), collect(forward_values)), "inverse-forward Laurent substitution")
            for variable in variables
        )...)
        return forward_then_inverse == variables && inverse_then_forward == variables
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _laurent_noether_endpoint_coefficient(entry, R, variable_index::Int, choose_degree)
    iszero(entry) && return zero(R)
    degrees = [Int(exponent_vector[variable_index]) for exponent_vector in collect(exponents(entry))]
    isempty(degrees) && return zero(R)
    endpoint_degree = choose_degree(degrees)
    total = zero(R)
    for (coefficient, exponent_vector) in zip(collect(coefficients(entry)), collect(exponents(entry)))
        exponent_vector[variable_index] == endpoint_degree || continue
        term = R(coefficient)
        for index in 1:ngens(R)
            index == variable_index && continue
            exponent = exponent_vector[index]
            exponent == 0 && continue
            term *= gen(R, index)^exponent
        end
        total += term
    end
    return total
end

_laurent_noether_leading_coefficient(entry, R, variable_index::Int) =
    _laurent_noether_endpoint_coefficient(entry, R, variable_index, maximum)

_laurent_noether_trailing_coefficient(entry, R, variable_index::Int) =
    _laurent_noether_endpoint_coefficient(entry, R, variable_index, minimum)

function _laurent_noether_is_supported_ring(R)::Bool
    _is_laurent_polynomial_ring(R) || return false
    ngens(R) == 2 || return false
    F = coefficient_ring(R)
    return Oscar.is_exact_type(typeof(zero(F))) && F isa Field
end

function _laurent_noether_validate_input(column::AbstractVector, selected_entry_index::Int, selected_generator)
    isempty(column) && throw(ArgumentError("Laurent Noether column must be nonempty"))
    R = _require_same_laurent_parent(column; label="Laurent Noether column")
    _laurent_noether_is_supported_ring(R) ||
        throw(ArgumentError("Laurent Noether certificate requires an exact field-backed Laurent ring with exactly two generators"))
    1 <= selected_entry_index <= length(column) ||
        throw(ArgumentError("selected Laurent column entry index is out of bounds"))
    variables = tuple(gens(R)...)
    selected_generator_index = findfirst(==(selected_generator), variables)
    selected_generator_index === nothing &&
        throw(ArgumentError("selected Laurent generator must be a generator of the Laurent ring"))
    iszero(column[selected_entry_index]) &&
        throw(ArgumentError("selected Laurent column entry must be nonzero"))
    return R, selected_generator_index, 3 - selected_generator_index
end

function _laurent_noether_power(entry, selected_generator_index::Int, other_generator_index::Int)::Int
    support = [Int.(collect(exponent_vector)) for exponent_vector in collect(exponents(entry))]
    isempty(support) && return 1
    selected_exponents = [exponent_vector[selected_generator_index] for exponent_vector in support]
    return max(1, maximum(selected_exponents) - minimum(selected_exponents) + 1)
end

function _laurent_noether_replay(column, R, selected_entry_index::Int, selected_generator_index::Int, other_generator_index::Int, noether_power::Int)
    variables = tuple(gens(R)...)
    forward_values = collect(variables)
    inverse_values = collect(variables)
    forward_values[other_generator_index] *= variables[selected_generator_index]^noether_power
    inverse_values[other_generator_index] *= variables[selected_generator_index]^(-noether_power)
    forward_substitution = _laurent_noether_map(variables, forward_values)
    inverse_substitution = _laurent_noether_map(variables, inverse_values)
    transformed_column = tuple((
        _coerce_into_ring(R, evaluate(entry, forward_values), "transformed Laurent column entry")
        for entry in column
    )...)
    replayed_selected_entry = transformed_column[selected_entry_index]
    leading_coefficient = _laurent_noether_leading_coefficient(replayed_selected_entry, R, selected_generator_index)
    trailing_coefficient = _laurent_noether_trailing_coefficient(replayed_selected_entry, R, selected_generator_index)
    return (; forward_substitution, inverse_substitution, transformed_column, replayed_selected_entry,
        leading_coefficient, trailing_coefficient,
        leading_coefficient_is_unit = is_unit(leading_coefficient),
        trailing_coefficient_is_unit = is_unit(trailing_coefficient))
end

function _laurent_noether_certificate(column::AbstractVector, selected_entry_index::Int, selected_generator)
    R, selected_generator_index, other_generator_index =
        _laurent_noether_validate_input(column, selected_entry_index, selected_generator)
    original_column = tuple(column...)
    noether_power = _laurent_noether_power(original_column[selected_entry_index], selected_generator_index, other_generator_index)
    replay = _laurent_noether_replay(original_column, R, selected_entry_index, selected_generator_index, other_generator_index, noether_power)
    certificate = LaurentNoetherCertificate(original_column, R, selected_entry_index,
        selected_generator_index, selected_generator, other_generator_index, noether_power,
        replay.forward_substitution, replay.inverse_substitution, replay.transformed_column,
        replay.replayed_selected_entry, replay.leading_coefficient, replay.trailing_coefficient,
        replay.leading_coefficient_is_unit, replay.trailing_coefficient_is_unit, :ok)
    status = _validate_laurent_noether_certificate(certificate)
    status == :ok || error("internal Laurent Noether certificate replay failed")
    return certificate
end

function _validate_laurent_noether_certificate(certificate)::Symbol
    try
        R = certificate.ring
        _laurent_noether_is_supported_ring(R) || return :invalid_ring
        column = certificate.original_column
        column isa Tuple || return :invalid_column
        isempty(column) && return :invalid_column
        all(entry -> parent(entry) === R, column) || return :wrong_ring
        1 <= certificate.selected_entry_index <= length(column) || return :invalid_selected_entry_index
        variables = tuple(gens(R)...)
        1 <= certificate.selected_generator_index <= 2 || return :invalid_selected_generator
        certificate.other_generator_index == 3 - certificate.selected_generator_index || return :invalid_generator_metadata
        certificate.selected_generator == variables[certificate.selected_generator_index] || return :invalid_generator_metadata
        expected_power = _laurent_noether_power(column[certificate.selected_entry_index], certificate.selected_generator_index, certificate.other_generator_index)
        certificate.noether_power == expected_power || return :invalid_noether_power
        replay = _laurent_noether_replay(column, R, certificate.selected_entry_index,
            certificate.selected_generator_index, certificate.other_generator_index, certificate.noether_power)
        _laurent_noether_maps_are_inverse(R, certificate.forward_substitution, certificate.inverse_substitution) || return :invalid_substitution_maps
        certificate.forward_substitution == replay.forward_substitution || return :invalid_substitution_maps
        certificate.inverse_substitution == replay.inverse_substitution || return :invalid_substitution_maps
        certificate.transformed_column == replay.transformed_column || return :stale_transformed_column
        certificate.replayed_selected_entry == replay.replayed_selected_entry || return :stale_selected_entry
        certificate.leading_coefficient == replay.leading_coefficient || return :stale_endpoint_coefficient
        certificate.trailing_coefficient == replay.trailing_coefficient || return :stale_endpoint_coefficient
        certificate.leading_coefficient_is_unit == replay.leading_coefficient_is_unit || return :stale_unit_status
        certificate.trailing_coefficient_is_unit == replay.trailing_coefficient_is_unit || return :stale_unit_status
        replay.leading_coefficient_is_unit && replay.trailing_coefficient_is_unit || return :nonunit_endpoint
        certificate.validation_status == :ok || return :invalid_validation_status
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_certificate
    end
end

struct LaurentToPolynomialColumnCertificate
    original_column
    ring
    noether_certificate
    selected_entry_index::Int
    selected_generator_index::Int
    selected_generator
    elementary_source_column
    conversion_factors::Vector
    conversion_product
    intermediate_laurent_column
    polynomial_ring
    polynomial_generators
    forward_polynomialization
    inverse_lift
    polynomial_column
    factor_lift_metadata
    replay
    validation_status::Symbol
end

struct LaurentToPolynomialECPBridgeCertificate
    conversion_certificate
    ordinary_child_certificate
    ordinary_factors::Vector
    raw_lifted_laurent_factors::Vector
    inverse_substituted_lifted_factors::Vector
    laurent_conversion_factors::Vector
    inverse_substituted_conversion_factors::Vector
    complete_factor_sequence::Vector
    target_basis_column
    recomputed_product
    replay_summary
    validation_status::Symbol
end

function _laurent_to_polynomial_validate_input(
    column::AbstractVector,
    noether_certificate,
    selected_entry_index::Int,
    selected_generator,
)
    length(column) >= 3 ||
        throw(ArgumentError("Laurent-to-polynomial conversion requires a column of length at least three"))
    R = _require_same_laurent_parent(column; label="Laurent-to-polynomial column")
    _laurent_noether_is_supported_ring(R) ||
        throw(ArgumentError("Laurent-to-polynomial conversion requires an exact field-backed Laurent ring with exactly two generators"))
    is_unimodular_column(column, R) ||
        throw(ArgumentError("Laurent-to-polynomial conversion requires a unimodular Laurent column"))
    _validate_laurent_noether_certificate(noether_certificate) == :ok ||
        throw(ArgumentError("Laurent-to-polynomial conversion requires a valid Laurent Noether certificate"))
    original_column = tuple(column...)
    noether_certificate.original_column == original_column ||
        throw(ArgumentError("Laurent Noether certificate is attached to a different source column"))
    noether_certificate.ring === R ||
        throw(ArgumentError("Laurent Noether certificate is attached to a different Laurent ring"))
    noether_certificate.selected_entry_index == selected_entry_index ||
        throw(ArgumentError("selected entry metadata does not match the Laurent Noether certificate"))
    variables = tuple(gens(R)...)
    selected_generator_index = findfirst(==(selected_generator), variables)
    selected_generator_index !== nothing ||
        throw(ArgumentError("selected Laurent generator must be a generator of the Laurent ring"))
    noether_certificate.selected_generator_index == selected_generator_index ||
        throw(ArgumentError("selected generator metadata does not match the Laurent Noether certificate"))
    noether_certificate.selected_generator == selected_generator ||
        throw(ArgumentError("selected generator does not match the Laurent Noether certificate"))
    _laurent_to_polynomial_is_supported_elementary_source_column(noether_certificate.transformed_column) ||
        throw(ArgumentError(
            "unsupported Laurent-to-polynomial conversion: Noether column must already be polynomial or contain a Laurent unit entry",
        ))
    return R, selected_generator_index
end

function _laurent_to_polynomial_column_matrix(column, R)
    return matrix(R, length(column), 1, collect(column))
end

function _laurent_to_polynomial_matrix_column_to_tuple(column_matrix)
    ncols(column_matrix) == 1 || throw(ArgumentError("expected a column matrix"))
    return tuple((column_matrix[row, 1] for row in 1:nrows(column_matrix))...)
end

function _laurent_to_polynomial_factor_sequence_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        _is_elementary_matrix_factor(factor, R, n) ||
            throw(ArgumentError("Laurent conversion factor must be elementary over the Laurent ring"))
        product *= factor
    end
    return product
end

function _laurent_to_polynomial_apply_factors(factors, column, R)
    product = _laurent_to_polynomial_factor_sequence_product(factors, R, length(column))
    return product * _laurent_to_polynomial_column_matrix(column, R)
end

function _laurent_to_polynomial_factor_sequences_equal(left, right)::Bool
    length(left) == length(right) || return false
    return all(index -> left[index] == right[index], eachindex(left))
end

function _laurent_to_polynomial_has_negative_exponent(entry)::Bool
    iszero(entry) && return false
    return any(exponent_vector -> any(<(0), Int.(collect(exponent_vector))), collect(exponents(entry)))
end

function _laurent_to_polynomial_is_polynomial_entry(entry)::Bool
    return !_laurent_to_polynomial_has_negative_exponent(entry)
end

function _laurent_to_polynomial_is_polynomial_column(column)::Bool
    return all(_laurent_to_polynomial_is_polynomial_entry, column)
end

function _laurent_to_polynomial_is_supported_elementary_source_column(column)::Bool
    return _laurent_to_polynomial_is_polynomial_column(column) || any(is_unit, column)
end

function _laurent_to_polynomial_unit_normalization_factors(n::Int, u, uinv, R)
    factors = typeof(identity_matrix(R, n))[]
    u == one(R) && return factors

    i = n - 1
    j = n
    for (row, col, coefficient) in (
        (i, j, u - one(R)),
        (j, i, one(R)),
        (i, j, uinv - one(R)),
        (j, i, -u),
    )
        push!(factors, elementary_matrix(n, row, col, coefficient, R))
    end
    return factors
end

function _laurent_to_polynomial_unit_preprocessing_factors(column, R)
    n = length(column)
    pivot_idx = findfirst(is_unit, column)
    pivot_idx === nothing &&
        throw(ArgumentError("unsupported Laurent-to-polynomial conversion: non-polynomial Noether column has no Laurent unit entry"))

    pivot = column[pivot_idx]
    pivot_inverse = inv(pivot)
    elimination_factors = typeof(identity_matrix(R, n))[]
    for row in 1:n
        row == pivot_idx && continue
        coefficient = -column[row] * pivot_inverse
        coefficient == zero(R) && continue
        push!(elimination_factors, elementary_matrix(n, row, pivot_idx, coefficient, R))
    end

    move_factors = typeof(identity_matrix(R, n))[]
    if pivot_idx != n
        push!(move_factors, elementary_matrix(n, pivot_idx, n, -one(R), R))
        push!(move_factors, elementary_matrix(n, n, pivot_idx, one(R), R))
    end

    normalization_factors = _laurent_to_polynomial_unit_normalization_factors(n, pivot, pivot_inverse, R)
    return vcat(normalization_factors, move_factors, elimination_factors)
end

function _laurent_to_polynomial_conversion_factors(elementary_source_column, R)
    factors = if _laurent_to_polynomial_is_polynomial_column(elementary_source_column)
        typeof(identity_matrix(R, length(elementary_source_column)))[]
    else
        _laurent_to_polynomial_unit_preprocessing_factors(elementary_source_column, R)
    end
    intermediate = _laurent_to_polynomial_matrix_column_to_tuple(
        _laurent_to_polynomial_apply_factors(factors, elementary_source_column, R),
    )
    _laurent_to_polynomial_is_polynomial_column(intermediate) ||
        throw(ArgumentError("unsupported Laurent-to-polynomial conversion: elementary preprocessing did not produce a polynomial column"))
    return factors
end

function _laurent_to_polynomial_target_ring(R)
    names = String.(string.(gens(R)))
    P, variables = suslin_polynomial_ring(coefficient_ring(R), names)
    return P, tuple(variables...)
end

function _laurent_to_polynomial_forward_map(R, P, polynomial_generators)
    laurent_generators = tuple(gens(R)...)
    return tuple((
        (; generator_index = index,
            laurent_generator = laurent_generators[index],
            polynomial_generator = polynomial_generators[index])
        for index in eachindex(laurent_generators)
    )...)
end

function _laurent_to_polynomial_inverse_lift_map(R, polynomial_generators)
    laurent_generators = tuple(gens(R)...)
    return tuple((
        (; generator_index = index,
            polynomial_generator = polynomial_generators[index],
            laurent_value = laurent_generators[index])
        for index in eachindex(laurent_generators)
    )...)
end

function _laurent_to_polynomial_factor_lift_metadata(R, P, inverse_lift, n::Int)
    return (;
        source_ring = R,
        polynomial_ring = P,
        inverse_lift,
        matrix_size = n,
        purpose = :ordinary_factor_lift_to_laurent_ring,
    )
end

function _laurent_to_polynomial_entry(entry, P, polynomial_generators)
    total = zero(P)
    iszero(entry) && return total
    for (coefficient, exponent_vector) in zip(collect(coefficients(entry)), collect(exponents(entry)))
        exponents_int = Int.(collect(exponent_vector))
        any(<(0), exponents_int) &&
            throw(ArgumentError("Laurent entry has a negative exponent and cannot be polynomialized"))
        term = P(coefficient)
        for index in eachindex(exponents_int)
            exponent = exponents_int[index]
            exponent == 0 && continue
            term *= polynomial_generators[index]^exponent
        end
        total += term
    end
    return total
end

function _laurent_to_polynomial_column(column, P, polynomial_generators)
    return tuple((_laurent_to_polynomial_entry(entry, P, polynomial_generators) for entry in column)...)
end

function _laurent_to_polynomial_lift_entry(entry, inverse_lift)
    isempty(inverse_lift) && throw(ArgumentError("inverse lift metadata must be nonempty"))
    R = parent(first(inverse_lift).laurent_value)
    total = zero(R)
    iszero(entry) && return total
    lift_values = tuple((map_entry.laurent_value for map_entry in inverse_lift)...)
    for (coefficient, exponent_vector) in zip(AbstractAlgebra.coefficients(entry), AbstractAlgebra.exponent_vectors(entry))
        term = R(coefficient)
        for index in eachindex(lift_values)
            exponent = Int(exponent_vector[index])
            exponent == 0 && continue
            term *= lift_values[index]^exponent
        end
        total += term
    end
    return total
end

function _laurent_to_polynomial_lift_factor(factor, factor_lift_metadata)
    P = factor_lift_metadata.polynomial_ring
    R = factor_lift_metadata.source_ring
    nrows(factor) == factor_lift_metadata.matrix_size ||
        throw(ArgumentError("ordinary factor has the wrong row count for Laurent lift metadata"))
    ncols(factor) == factor_lift_metadata.matrix_size ||
        throw(ArgumentError("ordinary factor has the wrong column count for Laurent lift metadata"))
    base_ring(factor) === P ||
        throw(ArgumentError("ordinary factor belongs to the wrong polynomial ring"))
    entries = [
        _laurent_to_polynomial_lift_entry(factor[row, col], factor_lift_metadata.inverse_lift)
        for col in 1:ncols(factor), row in 1:nrows(factor)
    ]
    return matrix(R, nrows(factor), ncols(factor), vec(entries))
end

function _laurent_to_polynomial_replay_summary(
    noether_status::Symbol,
    elementary_factors_are_elementary::Bool,
    elementary_replay_ok::Bool,
    polynomial_entries_in_ring::Bool,
    inverse_lift_ok::Bool,
    polynomial_unimodular::Bool,
)
    return (;
        noether_status,
        elementary_factors_are_elementary,
        elementary_replay_ok,
        polynomial_entries_in_ring,
        inverse_lift_ok,
        polynomial_unimodular,
        overall_ok = noether_status == :ok &&
            elementary_factors_are_elementary &&
            elementary_replay_ok &&
            polynomial_entries_in_ring &&
            inverse_lift_ok &&
            polynomial_unimodular,
    )
end

function _laurent_to_polynomial_expected_replay(certificate)
    R = certificate.ring
    P = certificate.polynomial_ring
    n = length(certificate.original_column)
    noether_status = _validate_laurent_noether_certificate(certificate.noether_certificate)
    elementary_factors_are_elementary = all(
        factor -> _is_elementary_matrix_factor(factor, R, n),
        certificate.conversion_factors,
    )
    replayed = _laurent_to_polynomial_matrix_column_to_tuple(
        _laurent_to_polynomial_apply_factors(
            certificate.conversion_factors,
            certificate.elementary_source_column,
            R,
        ),
    )
    elementary_replay_ok = replayed == certificate.intermediate_laurent_column
    polynomial_entries_in_ring = all(entry -> parent(entry) === P, certificate.polynomial_column)
    lifted_column = tuple((
        _laurent_to_polynomial_lift_entry(entry, certificate.inverse_lift)
        for entry in certificate.polynomial_column
    )...)
    inverse_lift_ok = lifted_column == certificate.intermediate_laurent_column
    polynomial_unimodular = polynomial_entries_in_ring &&
        is_unimodular_column(collect(certificate.polynomial_column), P)
    return _laurent_to_polynomial_replay_summary(
        noether_status,
        elementary_factors_are_elementary,
        elementary_replay_ok,
        polynomial_entries_in_ring,
        inverse_lift_ok,
        polynomial_unimodular,
    )
end

function _laurent_to_polynomial_certificate(
    column::AbstractVector,
    noether_certificate,
    selected_entry_index::Int,
    selected_generator,
)
    R, selected_generator_index = _laurent_to_polynomial_validate_input(
        column,
        noether_certificate,
        selected_entry_index,
        selected_generator,
    )
    original_column = tuple(column...)
    elementary_source_column = noether_certificate.transformed_column
    conversion_factors = _laurent_to_polynomial_conversion_factors(elementary_source_column, R)
    conversion_product =
        _laurent_to_polynomial_factor_sequence_product(conversion_factors, R, length(original_column))
    intermediate_laurent_column = _laurent_to_polynomial_matrix_column_to_tuple(
        conversion_product * _laurent_to_polynomial_column_matrix(elementary_source_column, R),
    )
    P, polynomial_generators = _laurent_to_polynomial_target_ring(R)
    forward_polynomialization =
        _laurent_to_polynomial_forward_map(R, P, polynomial_generators)
    inverse_lift = _laurent_to_polynomial_inverse_lift_map(R, polynomial_generators)
    polynomial_column =
        _laurent_to_polynomial_column(intermediate_laurent_column, P, polynomial_generators)
    factor_lift_metadata =
        _laurent_to_polynomial_factor_lift_metadata(R, P, inverse_lift, length(original_column))
    provisional = LaurentToPolynomialColumnCertificate(
        original_column,
        R,
        noether_certificate,
        selected_entry_index,
        selected_generator_index,
        selected_generator,
        elementary_source_column,
        conversion_factors,
        conversion_product,
        intermediate_laurent_column,
        P,
        polynomial_generators,
        forward_polynomialization,
        inverse_lift,
        polynomial_column,
        factor_lift_metadata,
        nothing,
        :ok,
    )
    replay = _laurent_to_polynomial_expected_replay(provisional)
    certificate = LaurentToPolynomialColumnCertificate(
        original_column,
        R,
        noether_certificate,
        selected_entry_index,
        selected_generator_index,
        selected_generator,
        elementary_source_column,
        conversion_factors,
        conversion_product,
        intermediate_laurent_column,
        P,
        polynomial_generators,
        forward_polynomialization,
        inverse_lift,
        polynomial_column,
        factor_lift_metadata,
        replay,
        :ok,
    )
    status = _validate_laurent_to_polynomial_certificate(certificate)
    status == :ok || error("internal Laurent-to-polynomial certificate replay failed with $(status)")
    return certificate
end

function _validate_laurent_to_polynomial_certificate(certificate)::Symbol
    try
        R = certificate.ring
        _laurent_noether_is_supported_ring(R) || return :invalid_ring
        original_column = certificate.original_column
        original_column isa Tuple || return :invalid_column
        length(original_column) >= 3 || return :invalid_column_length
        all(entry -> parent(entry) === R, original_column) || return :wrong_ring
        is_unimodular_column(collect(original_column), R) || return :nonunimodular_input

        noether_status = _validate_laurent_noether_certificate(certificate.noether_certificate)
        noether_status == :ok || return :invalid_noether_certificate
        certificate.noether_certificate.original_column == original_column ||
            return :invalid_noether_certificate
        certificate.noether_certificate.ring === R || return :invalid_noether_certificate

        variables = tuple(gens(R)...)
        1 <= certificate.selected_entry_index <= length(original_column) ||
            return :invalid_selected_entry_index
        1 <= certificate.selected_generator_index <= 2 ||
            return :invalid_selected_generator
        certificate.selected_generator == variables[certificate.selected_generator_index] ||
            return :invalid_generator_metadata
        certificate.noether_certificate.selected_entry_index == certificate.selected_entry_index ||
            return :invalid_noether_metadata
        certificate.noether_certificate.selected_generator_index == certificate.selected_generator_index ||
            return :invalid_noether_metadata
        certificate.noether_certificate.selected_generator == certificate.selected_generator ||
            return :invalid_noether_metadata

        elementary_source_column = certificate.noether_certificate.transformed_column
        certificate.elementary_source_column == elementary_source_column ||
            return :stale_elementary_source_column
        expected_factors = _laurent_to_polynomial_conversion_factors(elementary_source_column, R)
        _laurent_to_polynomial_factor_sequences_equal(certificate.conversion_factors, expected_factors) ||
            return :stale_conversion_factors
        n = length(original_column)
        expected_product =
            _laurent_to_polynomial_factor_sequence_product(certificate.conversion_factors, R, n)
        certificate.conversion_product == expected_product ||
            return :stale_conversion_product
        expected_intermediate = _laurent_to_polynomial_matrix_column_to_tuple(
            expected_product * _laurent_to_polynomial_column_matrix(elementary_source_column, R),
        )
        certificate.intermediate_laurent_column == expected_intermediate ||
            return :stale_intermediate_laurent_column

        P = certificate.polynomial_ring
        ngens(P) == ngens(R) || return :invalid_polynomial_ring
        coefficient_ring(P) == coefficient_ring(R) ||
            return :invalid_polynomial_ring
        polynomial_generators = tuple(gens(P)...)
        certificate.polynomial_generators == polynomial_generators ||
            return :invalid_polynomial_generators
        tuple(string.(polynomial_generators)...) == tuple(string.(gens(R))...) ||
            return :invalid_polynomial_ring

        expected_forward = _laurent_to_polynomial_forward_map(R, P, polynomial_generators)
        certificate.forward_polynomialization == expected_forward ||
            return :invalid_polynomialization_map
        expected_inverse = _laurent_to_polynomial_inverse_lift_map(R, polynomial_generators)
        certificate.inverse_lift == expected_inverse || return :invalid_inverse_lift
        expected_metadata =
            _laurent_to_polynomial_factor_lift_metadata(R, P, expected_inverse, n)
        certificate.factor_lift_metadata == expected_metadata ||
            return :invalid_factor_lift_metadata

        expected_polynomial_column =
            _laurent_to_polynomial_column(expected_intermediate, P, polynomial_generators)
        certificate.polynomial_column == expected_polynomial_column ||
            return :stale_polynomial_column
        all(entry -> parent(entry) === P, certificate.polynomial_column) ||
            return :wrong_polynomial_ring
        lifted_column = tuple((
            _laurent_to_polynomial_lift_entry(entry, certificate.inverse_lift)
            for entry in certificate.polynomial_column
        )...)
        lifted_column == expected_intermediate || return :invalid_inverse_lift
        is_unimodular_column(collect(certificate.polynomial_column), P) ||
            return :polynomial_column_not_unimodular

        expected_replay = _laurent_to_polynomial_expected_replay(certificate)
        certificate.replay == expected_replay || return :stale_replay
        expected_replay.overall_ok || return :invalid_replay
        certificate.validation_status == :ok || return :invalid_validation_status
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_certificate
    end
end

function _laurent_to_polynomial_target_basis_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _laurent_to_polynomial_ecp_bridge_product(factors, R, n::Int)
    return _laurent_to_polynomial_factor_sequence_product(factors, R, n)
end

function _laurent_to_polynomial_ecp_bridge_lifted_factors(conversion_certificate, ordinary_factors)
    return [
        _laurent_to_polynomial_lift_factor(factor, conversion_certificate.factor_lift_metadata)
        for factor in ordinary_factors
    ]
end

function _laurent_to_polynomial_ecp_bridge_inverse_substituted_factors(factors, conversion_certificate)
    return _ecp_substitute_factor_sequence(
        factors,
        conversion_certificate.noether_certificate.inverse_substitution,
        conversion_certificate.ring,
    )
end

function _laurent_to_polynomial_ecp_bridge_replay_summary(certificate)
    conversion_certificate = certificate.conversion_certificate
    conversion_status = _validate_laurent_to_polynomial_certificate(conversion_certificate)
    R = conversion_certificate.ring
    n = length(conversion_certificate.original_column)
    target = _laurent_to_polynomial_target_basis_column(R, n)

    ordinary_child_ok = verify_ecp_column_reduction(certificate.ordinary_child_certificate)
    ordinary_child_column_ok = ordinary_child_ok &&
        certificate.ordinary_child_certificate.original_column ==
        collect(conversion_certificate.polynomial_column)
    ordinary_child_ring_ok = ordinary_child_ok &&
        certificate.ordinary_child_certificate.ring === conversion_certificate.polynomial_ring
    ordinary_factors_ok = ordinary_child_ok &&
        certificate.ordinary_factors == certificate.ordinary_child_certificate.factors

    expected_raw_lifted = ordinary_factors_ok ?
        _laurent_to_polynomial_ecp_bridge_lifted_factors(
            conversion_certificate,
            certificate.ordinary_child_certificate.factors,
        ) :
        Any[]
    raw_lifted_factors_ok =
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.raw_lifted_laurent_factors,
            expected_raw_lifted,
        )
    raw_lifted_factors_are_elementary = all(
        factor -> _is_elementary_matrix_factor(factor, R, n),
        certificate.raw_lifted_laurent_factors,
    )

    expected_inverse_lifted = raw_lifted_factors_ok ?
        _laurent_to_polynomial_ecp_bridge_inverse_substituted_factors(
            expected_raw_lifted,
            conversion_certificate,
        ) :
        Any[]
    inverse_substituted_lifted_factors_ok =
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.inverse_substituted_lifted_factors,
            expected_inverse_lifted,
        )
    inverse_substituted_lifted_factors_are_elementary = all(
        factor -> _is_elementary_matrix_factor(factor, R, n),
        certificate.inverse_substituted_lifted_factors,
    )

    laurent_conversion_factors_ok =
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.laurent_conversion_factors,
            conversion_certificate.conversion_factors,
        )
    expected_inverse_conversion = laurent_conversion_factors_ok ?
        _laurent_to_polynomial_ecp_bridge_inverse_substituted_factors(
            conversion_certificate.conversion_factors,
            conversion_certificate,
        ) :
        Any[]
    inverse_substituted_conversion_factors_ok =
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.inverse_substituted_conversion_factors,
            expected_inverse_conversion,
        )
    inverse_substituted_conversion_factors_are_elementary = all(
        factor -> _is_elementary_matrix_factor(factor, R, n),
        certificate.inverse_substituted_conversion_factors,
    )

    expected_complete_sequence = vcat(expected_inverse_lifted, expected_inverse_conversion)
    complete_factor_sequence_ok =
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.complete_factor_sequence,
            expected_complete_sequence,
        )
    complete_factors_are_elementary = all(
        factor -> _is_elementary_matrix_factor(factor, R, n),
        certificate.complete_factor_sequence,
    )

    expected_product = complete_factor_sequence_ok ?
        _laurent_to_polynomial_ecp_bridge_product(expected_complete_sequence, R, n) :
        identity_matrix(R, n)
    recomputed_product_ok = certificate.recomputed_product == expected_product
    target_basis_column_ok = certificate.target_basis_column == target
    replayed_column = certificate.recomputed_product *
        _laurent_to_polynomial_column_matrix(conversion_certificate.original_column, R)
    replay_ok = replayed_column == target

    overall_ok = conversion_status == :ok &&
        ordinary_child_ok &&
        ordinary_child_column_ok &&
        ordinary_child_ring_ok &&
        ordinary_factors_ok &&
        raw_lifted_factors_ok &&
        raw_lifted_factors_are_elementary &&
        inverse_substituted_lifted_factors_ok &&
        inverse_substituted_lifted_factors_are_elementary &&
        laurent_conversion_factors_ok &&
        inverse_substituted_conversion_factors_ok &&
        inverse_substituted_conversion_factors_are_elementary &&
        complete_factor_sequence_ok &&
        complete_factors_are_elementary &&
        recomputed_product_ok &&
        target_basis_column_ok &&
        replay_ok

    return (;
        conversion_status,
        ordinary_child_ok,
        ordinary_child_column_ok,
        ordinary_child_ring_ok,
        ordinary_factors_ok,
        raw_lifted_factors_ok,
        raw_lifted_factors_are_elementary,
        inverse_substituted_lifted_factors_ok,
        inverse_substituted_lifted_factors_are_elementary,
        laurent_conversion_factors_ok,
        inverse_substituted_conversion_factors_ok,
        inverse_substituted_conversion_factors_are_elementary,
        complete_factor_sequence_ok,
        complete_factors_are_elementary,
        recomputed_product_ok,
        target_basis_column_ok,
        replay_ok,
        overall_ok,
        replayed_column,
    )
end

function _laurent_to_polynomial_ecp_bridge_child_certificate(conversion_certificate)
    child_error = nothing
    ordinary_child_certificate = try
        ecp_column_reduction_certificate(
            collect(conversion_certificate.polynomial_column),
            conversion_certificate.polynomial_ring,
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        child_error = err
        nothing
    end
    ordinary_child_certificate !== nothing && return ordinary_child_certificate

    fallback_child_certificate = _ecp_rank_one_normality_unit_child_certificate(
        collect(conversion_certificate.polynomial_column),
        conversion_certificate.polynomial_ring,
    )
    fallback_child_certificate !== nothing && return fallback_child_certificate
    throw(child_error)
end

function _laurent_to_polynomial_ecp_bridge_certificate(conversion_certificate)
    status = _validate_laurent_to_polynomial_certificate(conversion_certificate)
    status == :ok ||
        throw(ArgumentError("Laurent-to-polynomial ECP bridge requires a verified conversion certificate"))

    ordinary_child_certificate =
        _laurent_to_polynomial_ecp_bridge_child_certificate(conversion_certificate)
    verify_ecp_column_reduction(ordinary_child_certificate) ||
        throw(ArgumentError("ordinary ECP child certificate does not verify"))

    ordinary_factors = collect(ordinary_child_certificate.factors)
    raw_lifted_laurent_factors =
        _laurent_to_polynomial_ecp_bridge_lifted_factors(
            conversion_certificate,
            ordinary_factors,
        )
    inverse_substituted_lifted_factors =
        _laurent_to_polynomial_ecp_bridge_inverse_substituted_factors(
            raw_lifted_laurent_factors,
            conversion_certificate,
        )
    laurent_conversion_factors = collect(conversion_certificate.conversion_factors)
    inverse_substituted_conversion_factors =
        _laurent_to_polynomial_ecp_bridge_inverse_substituted_factors(
            laurent_conversion_factors,
            conversion_certificate,
        )
    complete_factor_sequence = vcat(
        inverse_substituted_lifted_factors,
        inverse_substituted_conversion_factors,
    )
    R = conversion_certificate.ring
    n = length(conversion_certificate.original_column)
    target_basis_column = _laurent_to_polynomial_target_basis_column(R, n)
    recomputed_product =
        _laurent_to_polynomial_ecp_bridge_product(complete_factor_sequence, R, n)
    provisional = LaurentToPolynomialECPBridgeCertificate(
        conversion_certificate,
        ordinary_child_certificate,
        ordinary_factors,
        raw_lifted_laurent_factors,
        inverse_substituted_lifted_factors,
        laurent_conversion_factors,
        inverse_substituted_conversion_factors,
        complete_factor_sequence,
        target_basis_column,
        recomputed_product,
        nothing,
        :ok,
    )
    replay_summary = _laurent_to_polynomial_ecp_bridge_replay_summary(provisional)
    replay_summary.overall_ok ||
        error("internal Laurent-to-polynomial ECP bridge replay failed")
    certificate = LaurentToPolynomialECPBridgeCertificate(
        conversion_certificate,
        ordinary_child_certificate,
        ordinary_factors,
        raw_lifted_laurent_factors,
        inverse_substituted_lifted_factors,
        laurent_conversion_factors,
        inverse_substituted_conversion_factors,
        complete_factor_sequence,
        target_basis_column,
        recomputed_product,
        replay_summary,
        :ok,
    )
    _validate_laurent_to_polynomial_ecp_bridge_certificate(certificate) == :ok ||
        error("internal Laurent-to-polynomial ECP bridge certificate verification failed")
    return certificate
end

function _validate_laurent_to_polynomial_ecp_bridge_certificate(certificate)::Symbol
    try
        conversion_status =
            _validate_laurent_to_polynomial_certificate(certificate.conversion_certificate)
        conversion_status == :ok || return :invalid_conversion_certificate
        verify_ecp_column_reduction(certificate.ordinary_child_certificate) ||
            return :invalid_ordinary_child_certificate
        certificate.ordinary_child_certificate.original_column ==
            collect(certificate.conversion_certificate.polynomial_column) ||
            return :stale_ordinary_child_certificate
        certificate.ordinary_child_certificate.ring ===
            certificate.conversion_certificate.polynomial_ring ||
            return :stale_ordinary_child_certificate
        certificate.ordinary_factors == certificate.ordinary_child_certificate.factors ||
            return :stale_ordinary_factors

        expected_raw_lifted =
            _laurent_to_polynomial_ecp_bridge_lifted_factors(
                certificate.conversion_certificate,
                certificate.ordinary_child_certificate.factors,
            )
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.raw_lifted_laurent_factors,
            expected_raw_lifted,
        ) || return :stale_lifted_laurent_factors

        expected_inverse_lifted =
            _laurent_to_polynomial_ecp_bridge_inverse_substituted_factors(
                expected_raw_lifted,
                certificate.conversion_certificate,
            )
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.inverse_substituted_lifted_factors,
            expected_inverse_lifted,
        ) || return :stale_inverse_substituted_lifted_factors

        _laurent_to_polynomial_factor_sequences_equal(
            certificate.laurent_conversion_factors,
            certificate.conversion_certificate.conversion_factors,
        ) || return :stale_laurent_conversion_factors
        expected_inverse_conversion =
            _laurent_to_polynomial_ecp_bridge_inverse_substituted_factors(
                certificate.conversion_certificate.conversion_factors,
                certificate.conversion_certificate,
            )
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.inverse_substituted_conversion_factors,
            expected_inverse_conversion,
        ) || return :stale_inverse_substituted_conversion_factors

        expected_complete_sequence = vcat(expected_inverse_lifted, expected_inverse_conversion)
        _laurent_to_polynomial_factor_sequences_equal(
            certificate.complete_factor_sequence,
            expected_complete_sequence,
        ) || return :stale_complete_factor_sequence

        R = certificate.conversion_certificate.ring
        n = length(certificate.conversion_certificate.original_column)
        expected_target = _laurent_to_polynomial_target_basis_column(R, n)
        certificate.target_basis_column == expected_target ||
            return :stale_target_basis_column
        expected_product =
            _laurent_to_polynomial_ecp_bridge_product(expected_complete_sequence, R, n)
        certificate.recomputed_product == expected_product ||
            return :stale_recomputed_product

        expected_replay = _laurent_to_polynomial_ecp_bridge_replay_summary(certificate)
        expected_replay.overall_ok || return :invalid_replay
        certificate.replay_summary == expected_replay || return :stale_replay
        certificate.validation_status == :ok || return :invalid_validation_status
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_certificate
    end
end
