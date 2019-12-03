include("lanczos.jl")

using LinearAlgebra

"""
Lanczos-based MINRES.
"""
function minres(A, b, maxit)
    n = length(b)
    x = zeros(n)
    Q = zeros(n, maxit + 2)
    H = zeros(maxit + 2, maxit + 1)

    β = norm(b)
    r = β # will become norm(b - Ax_k)
    # finish

    Q[:, 2] = b / norm(b)
    for k = 2 : maxit + 1
        _q = A * Q[:, k]
        q = _q - Q[:, k]' * _q * Q[:, k] - Q[:, k-1]' * _q * Q[:, k-1]
        H[k + 1, k] = norm(q)
        H[k, k] = Q[:, k]' * _q
        H[k - 1, k] = Q[:, k - 1] * _q
        Q[:, k + 1] = q / norm(q)


    end

    return x, r
end
