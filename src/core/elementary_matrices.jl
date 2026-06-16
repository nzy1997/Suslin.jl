function elementary_matrix(n::Int, i::Int, j::Int, a, R)
    i == j && throw(ArgumentError("i and j must differ"))

    E = identity_matrix(R, n)
    coerced_a = try
        R(a)
    catch
        throw(ArgumentError("a must be coercible into the target ring"))
    end
    E[i, j] = coerced_a
    return E
end
