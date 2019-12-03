include("lanczos.jl")

using LinearAlgebra

"""
Simple implementation of the Conjugate Gradient Method for a symmetric
positive definite matrix A.
"""
function cg(A, b, k)
    V, T, ρ = lanczos(A, b, k)
    rhs = zeros(k)
    rhs[1] = norm(b)
    y = T \ rhs
    return V * y # appromixation for x
end
