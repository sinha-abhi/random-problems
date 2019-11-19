include("lanczos.jl")

using LinearAlgebra
using Plots

"""
Generates an n x n Strakos matrix.

A Strakos matrix is a diagonal matrix, with the ith element being
    d_i = λ_1 + (i - 1) / (n -1) (λ_n - λ_1) ρ^(n - i).
"""
function strakos(λ_1, λ_n, ρ, n)
    δλ = λ_n - λ_1
    d = [λ_1 + (i - 1) / (n - 1) * δλ * ρ^(n - i) for i in collect(1 : n)]
    return Diagonal(d)
end

function residual(V, matrix=true)
    return log(norm(V' * V - I) + 10^(-20))
end

function col_residual(V, offset=false)
    k = size(V, 2)
    if !offset
        v1 = V[:, 1]
        return [log(abs(v1' * V[:, i]) + 10^(-20)) for i in collect(1 : k-1)]
    else
        @assert(k > 2, "must have at least 3 columns in V")
        return [log(abs(V[:, i-2]' * V[:, i]) + 10^(-20)) for i in collect(3 : k-1)]
    end
end

##
"""
Lanczos method on Strakos matrices.
"""
n =  30
v1 = ones(n) / sqrt(n) # v1 = e / √n
rv = rand(n)


R = zeros(n, 2) # plot info
S = strakos(0.1, 100, 0.9, n)
for k = 1 : n
    V, _, _ = lanczos(S, v1, k)
    U, _, _ = lanczos(S, rv, k)
    # matrix residuals
    R[k, 1] = residual(V)
    R[k, 2] = residual(U)
end

plot(1 : n, R, label = ["v1" "rv"], title = "Matrix Residuals")

V, _, _ = lanczos(S, v1, n)
U, _, _ = lanczos(S, rv, n)
CR = [col_residual(V) col_residual(U)] # column residuals
OCR = [col_residual(V, true) col_residual(U, true)] # offset column residuals
plot(1 : n, CR, label = ["v1" "rv"], title = "Column Residuals")
plot(1 : n - 2, OCR, label = ["v1" "rv"], title = "Offset Column residuals")
