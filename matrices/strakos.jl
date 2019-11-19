include("lanczos.jl")

using LinearAlgebra
using Plots

"""
Generates an n x n Strakos matrix.

A Strakos matrix is a diagonal matrix, with the ith element being
    d_i = λ_1 + (i - 1) / (n -1) (λ_n - λ_1) ρ^(n - i).
"""
function strakos(λ_1, λ_n, ρ, n)
    d = zeros(n)
    δλ = λ_n - λ_1
    for i = 1 : n
        d[i] = λ_1 + (i - 1) / (n - 1) * δλ * ρ^(n - i)
    end

    return Diagonal(d)
end

function residual(V)
    return log(norm(V' * V - I) + 10^(-20))
end

## Lanczos method on Strakos matrices
n =  30
v1 = ones(n) / sqrt(n) # v1 = e / √n
rv = rand(n)

R = zeros(n, 2)
S = strakos(0.1, 100, 0.9, n)
for k = 1 : n
    T, ρ = lanczos(S, v1, k)
    R[k, 1] = residual(T)
    U, ξ = lanczos(S, rv, k)
    R[k, 2] = residual(U)
end

plot(1 : n, R, label = ["v1" "rv"])
