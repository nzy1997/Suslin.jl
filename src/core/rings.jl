function suslin_polynomial_ring(F, names::Vector{String})
    R, vars = Oscar.polynomial_ring(F, names)
    return R, collect(vars)
end
