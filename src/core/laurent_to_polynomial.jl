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
