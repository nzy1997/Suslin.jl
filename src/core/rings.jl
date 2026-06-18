function suslin_polynomial_ring(F, names::Vector{String})
    R, vars = Oscar.polynomial_ring(F, names)
    return R, collect(vars)
end

function suslin_laurent_polynomial_ring(F, names::Vector{String})
    R, vars = Oscar.laurent_polynomial_ring(F, names)
    return R, collect(vars)
end

function _is_laurent_polynomial_ring(R)
    return R isa LaurentPolyRing || R isa LaurentMPolyRing
end

function _require_laurent_polynomial_ring(R; label::AbstractString="ring")
    _is_laurent_polynomial_ring(R) && return R
    throw(ArgumentError("$label must be a Laurent polynomial ring"))
end

function _parent_for_validation(value, label::AbstractString)
    try
        return parent(value)
    catch err
        err isa MethodError || rethrow()
        throw(ArgumentError("$label must be an element of a Laurent polynomial ring"))
    end
end

function _require_laurent_element(value; label::AbstractString="value")
    _require_laurent_polynomial_ring(_parent_for_validation(value, label); label="$label parent")
    return value
end

function _require_laurent_element(value, R; label::AbstractString="value")
    _require_laurent_polynomial_ring(R; label="expected parent")
    parent_value = _parent_for_validation(value, label)
    parent_value === R || throw(ArgumentError("$label must belong to the expected Laurent polynomial ring"))
    return value
end

function _require_same_laurent_parent(values; label::AbstractString="values")
    state = iterate(values)
    state === nothing && throw(ArgumentError("$label must be nonempty"))

    first_value, next_state = state
    R = _require_laurent_polynomial_ring(_parent_for_validation(first_value, label); label="$label parent")

    while true
        state = iterate(values, next_state)
        state === nothing && return R
        value, next_state = state
        _require_laurent_element(value, R; label=label)
    end
end
