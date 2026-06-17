struct LocalCertificate
    indices::Vector{Int}
    denominators::Vector

    function LocalCertificate(indices::Vector{Int}, denominators::Vector)
        length(indices) == length(denominators) || throw(ArgumentError("indices and denominators must have the same length"))
        return new(indices, denominators)
    end
end

function common_denominator_factor(entries::AbstractVector)
    Base.require_one_based_indexing(entries)
    isempty(entries) && throw(ArgumentError("entries must not be empty"))

    factor = _denominator_factor(entries[1])
    for idx in 2:length(entries)
        factor *= _denominator_factor(entries[idx])
    end
    return factor
end

function _denominator_factor(entry)
    try
        return denominator(entry)
    catch err
        if err isa MethodError
            return one(parent(entry))
        end
        rethrow()
    end
end
