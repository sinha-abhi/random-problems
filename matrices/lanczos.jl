using LinearAlgebra

"""
Implementation of the Lanczos method to compute the orthonormal basis for the
Krylov space of a Hermitian matrix.
"""
function lanczos(A::AbstractMatrix, b::AbstractVector, k)
    n = size(A, 1)
    V = zeros(n, k + 1)
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
