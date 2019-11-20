using LinearAlgebra

"""
Lanczos method to compute the orthonormal basis for the Krylov space
of a symmetric matrix.

Julia doesn't allow for non-square Tridiagonal matrices, so ρ takes the place
of the missing element of the last row.
"""
function lanczos(A::AbstractMatrix, b::AbstractVector, k)
    n = size(A, 1)
    V = zeros(n, k + 1)
    T = Tridiagonal(zeros(k - 1), zeros(k), zeros(k - 1))
    ρ = 0.0
    V[:, 1] = b / norm(b)
    for i = 1 : k
        y = A * V[:, i]
        for j = max(i - 1, 1) : i
            T[j, i] = V[:, j]' * y
            y -= T[j, i] * V[:, i]
        end
        ρ = norm(y)
        V[:, i + 1] = y / ρ
        if i < k
            T[i + 1, i] = ρ
        end
    end

    return V, T, ρ
end

"""
More efficient implementation of Lanczos.
"""
function eff_lanczos(A::AbstractMatrix, b::AbstractVector, k)
    n = size(A, 1)
    T = Tridiagonal(zeros(k - 1), zeros(k), zeros(k - 1))
    ρ = 0.0
    vold = 0.0
    v = b / norm(b)
    for j = 1 : k
        y = A * v
        if j > 1
            T[j - 1, j] = vold' * y
            y -= T[j - 1, j] * vold
        end
        T[j, j] = v' * y
        y -= T[j, j] * v

        vold = v
        ρ = norm(y)
        v = y / ρ
        if j < k
            T[j + 1, j] = ρ
        end
    end

    return T, ρ
end

"""
Computes the next iteration of the Lanczos method.

Let u be v_(k-1), v be v_k, and β be β_k. We return α_k, and β, v for the
(k+1)st iteration.
"""
function lanczos_iter(A, u, v, β)
    y = A * v
    α = v' * y # compute α_k
    y -= α * v
    _v = y - β * u # compute v_(k+1)
    _β = norm(_v)
    if _β != 0
        _v /= _β
    end

    return α, _β, _v
end
