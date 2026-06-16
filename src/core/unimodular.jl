function is_unimodular_column(v::AbstractVector, R)
    Base.require_one_based_indexing(v)
    isempty(v) && throw(ArgumentError("v must be nonempty"))

    column = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in eachindex(v)]
    return contains(ideal(R, column), ideal(R, [one(R)]))
end
