using LinearAlgebra

"""
Searches the Krylov subspace of a symmetric matrix A of size nxn
the fact that λ_max = max{x^TAx}, where x is a unit eigenvector.
"""
function eigen_max(A::AbstractMatrix, b::AbstractVector)
    @assert(issymmetric(A), "A must be symmetric")
    @assert(size(A, 1) == size(b, 1), "vector b has incorrect dimensions")

    ε = 1e-12
    Q, H = arnoldi(A, b, size(A, 1))
    Λ = zeros(size(Q, 2))

    for i = 1 : size(Q, 2)
        q = Q[:, i]
        if norm(q - zeros(size(q))) < ε
            continue
        end
        α = q' * A * q
        println(α)
        Λ[i] = α
    end

    return maximum(Λ)
end

function arnoldi(A::AbstractMatrix, b::AbstractVector, k)
    n = size(A, 1)
    V = zeros(n, k + 1) # orthonormal basis for Krylov subspace
    H = zeros(k + 1, k) # basis for A on V

    Q[:, 1] = b / norm(b)
    for j = 1 : k
        y = A * V[:, j]
        for i = 1 : k
            H[i, j] = V[:, i]' * y
            y -= H[i, j] * V[:, i]
        end
        H[j + 1, j] = norm(y)
        V[:, j + 1] = y / H[j + 1, j]
    end

    return V, H
end

## An example test case
A = [[1 2 3]; [2 3 4]; [3 4 5]]
b = [1; 1; 1]
Q, H = arnoldi(A, b, size(A, 1))
λ = eigen_max(A, b)
α = eigmax(A)
println(λ)
print(α)
@show norm(α - λ)
